clear all
close all
clc

%% Variables and constants
mpump = 50; 
Volume = 84;
DensiteAir = 1.204;
e = 0.355;
k = 2781;
s1 = (6*4) * 2;
s2 = (6*3.5) * 2;
s3 = (4*3.5) * 2;
s = s1 + s2 + s3;
R = e * k/s ;
M = Volume * DensiteAir;
c = 1256;
sfe = 1.5;
Q = 100;
Az  = (-(1/M)*mpump)-1/(M*c*R);
theta = 0.8;
a = 1;
beta = 1;
sigma = 0.8;
k = 1000;

%% Read data from files
profile = dlmread('profile.txt');
data = dlmread('data.txt');
data(1, :) = [];
%% Mpc Parameters
ts = 0.2;
Tstop = 72;
p = 24;
m = 3;
time = (0 : ts : (Tstop)-ts);

%% Jacobians matrices definitions
A = [Az (mpump/M) 0;
       theta a 0;
       0 0 1];
   
B = [0 0 0 (1/M*c*R) (sfe/M*c) (Q/M*c);
       beta 0 0 sigma 0 0;
       0 ts 0 0 0 0];

C = [1 0 0;
       0 0 1;
       0 0 0];

D = [0 0 0 0 0 0;
       0 1 0 0 0 0;
       k ts ts 0 0 0];
   
plant = ss(A, B, C, D);
plant.InputDelay = 0.001;

%% Reference definition
Tzr = [zeros(1, 35) ones(1, 25) zeros(1, 5) ones(1, 30) zeros(1, 30)];
Tzr = [Tzr Tzr Tzr 0];
Tzr = 21 *Tzr;
Tzr(length(time)+1:end)= [];

%% Measurables disturbances
Ppv = data(:, 3);% Ppv : Production panneaux photovoltaiques pendant 3 jours en Kj/h
I = data(:, 2); % Radiation solaire pendant 3 jours 
Tout = data(:, 4); % temperature ambiante
x = 0:71;
Tout = interp1(x', Tout, time);
I = interp1(x', I, time);
Ppv = interp1(x', Ppv, time);

for i = 1 : length(Tout)
    if isnan(Ppv(i))
        Ppv(i) = Ppv(i-1);
    end
    if isnan(Tout(i))
        Tout(i) = Tout(i-1);
    end
    if isnan(I(i))
        I(i) = I(i-1);
    end
end

v = [Tout' zeros(length(time), 2)];

%% MPC creation
plant = setmpcsignals(plant, 'MV', 1:3, 'MD', 4:6);
mpcBc = mpc(plant, ts, p, m);
% mpcBc.Model.Nominal.X= [15, 0, 0];
% pl = [320*ones(1, 9),  250, 400*ones(1, 10), 360*ones(1, 4)];
r = [Tzr; zeros(1, 360) ;Ppv];
% 
%% Scale 
% manipulated variables
mpcBc.ManipulatedVariables(1).ScaleFactor = 40;
mpcBc.ManipulatedVariables(2).ScaleFactor = 1000;
mpcBc.ManipulatedVariables(3).ScaleFactor = 1000;
% disturbance 
mpcBc.DisturbanceVariables(1).ScaleFactor = 40;
% Output
mpcBc.OutputVariables(1).ScaleFactor = 40;
mpcBc.OutputVariables(2).ScaleFactor = 1000;
mpcBc.OutputVariables(3).ScaleFactor = 1000;

%% Weigths
% Manipulated variables
mpcBc.Weights.ManipulatedVariables = [1 1 1];
% Output Variables 
mpcBc.Weights.OutputVariables = [1 0 1];

% %% constraints
% % Manipulated variables
% mpcBc.ManipulatedVariables(1).Max = 60;
% mpcBc.ManipulatedVariables(1).RateMax = 15;
% mpcBc.ManipulatedVariables(1).Min = 0;
% mpcBc.ManipulatedVariables(1).RateMin = -15;
% mpcBc.ManipulatedVariables(1).MaxECR = 3;
% mpcBc.ManipulatedVariables(1).MinECR = 1;
% 
% mpcBc.ManipulatedVariables(2).Max = 8000;
% mpcBc.ManipulatedVariables(2).RateMax = 1000;
% mpcBc.ManipulatedVariables(2).Min = -8000;
% mpcBc.ManipulatedVariables(2).RateMin = -1000;
% mpcBc.ManipulatedVariables(2).MaxECR = 3;
% mpcBc.ManipulatedVariables(2).MinECR = 1;
% 
% mpcBc.ManipulatedVariables(3).Max = 8000;
% mpcBc.ManipulatedVariables(3).RateMax = 1000;
% mpcBc.ManipulatedVariables(3).Min = -8000;
% mpcBc.ManipulatedVariables(3).RateMin = -1000;
% mpcBc.ManipulatedVariables(3).MaxECR = 3;
% mpcBc.ManipulatedVariables(3).MinECR = 1;
% 
% % outputs variables
% mpcBc.OutputVariables(1).Max = 25;
% mpcBc.OutputVariables(1).MaxECR = 0;
% mpcBc.OutputVariables(1).Min = -25;
% mpcBc.OutputVariables(1).MinECR = 0;
% 
% mpcBc.OutputVariables(2).Max = 10000;
% mpcBc.OutputVariables(2).MaxECR = 3;
% mpcBc.OutputVariables(2).Min = 1000;
% mpcBc.OutputVariables(2).MinECR = 1;
% 
% 
% params = mpcsimopt(mpcBc);
% params.MDLookAhead = 'on';
% params.RefLookAhead = 'on';
% [y, t, u] = sim(mpcBc, Tstop/ts, r',v, params);
% sim(mpcBc, Tstop/ts, r',v, params);
% 
% 
