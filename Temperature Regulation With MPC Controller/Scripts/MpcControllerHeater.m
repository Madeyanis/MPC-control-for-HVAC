clear all
clc
close all


%% Cette experience est realisee pour une duree de 3 jours [1 janvier - 4 janvier[


%% Problem definition 
% dx/dt = Ax + Bu + Bd
% y = Cx
% x = [Tz] state représente la température de la zone
% y = [Tz] sortie     //     //    //       // //  //
% u = [Theat] control signal: température du chauffage
% d = [TOut Ie Iw Is In N] les entrées non commandables mais mesurables

%% Load fichier
data = readtable('data.txt'); data = table2array(data);
data(1, :) = [];
%% constants definition and load data
Volume = 84;
DensiteAir = 1.204;
% M = Masse d'air à l'interieur de la zone (kg)
M = Volume * DensiteAir;
% c = Capacité calorifique de l'air à pression constante (j/kg.K)
c = 1256; c = c / 3600;
% m = Débit du fluid entrant dans le heater
mpump = 50;
% Theat = température du heater (commande à réguler)
% Q = Apport d'un seul occupant (watt)
Q = 150;
% N = Nombre d'occupant
N = 2 * [zeros(1, 39), ones(1, 21), zeros(1, 5), ones(1, 25), zeros(1, 30)];
N = repmat(N, 1, 3);
% Sf = surface exposé au soleil
Sf = 1.5;
% R = résistance thermique de l'enveloppe du bâtiment, R = e / (k*s)
% k = conductivité thermique du matériau, S = surface du materiau
% e = epaisseur du matériau
e = 0.355;
k = 2781;
s1 = (6*4) * 2;
s2 = (6*3.5) * 2;
s3 = (4*3.5) * 2;
s = s1 + s2 + s3;
R = e * k/s ;
R = R*1;
% Ppv : Production panneaux photovoltaiques pendant 3 jours
Ppv = data(:, 2); Ppv = Ppv * 1000;
% Radiation solaire pendant 3 jours 
I = data(:, 4); I = I * 1000;
% temperature ambiante
AmbianTemperature = data(:, 5);
% Chauffage
cpf = (4.19/3.6);

%% Some variables for Mpc
Ts = 0.2;
P = 8;
m = 3;
Tstop = 72;

time = (0 : Ts : (Tstop)-Ts);

% Reference definition
r = [zeros(1, 35) ones(1, 25) zeros(1, 5) ones(1, 30) zeros(1, 30)];
r = [r r r 0];
r = 23 *r;
r(length(time)+1:end)= [];

x = 0:71;
AmbianTemperatureInterpol = interp1(x', AmbianTemperature, time);

for i = 1 : length(AmbianTemperatureInterpol)

    if isnan(AmbianTemperatureInterpol(i))
        AmbianTemperatureInterpol(i) = AmbianTemperatureInterpol(i-1);
    end
    
end

I2 = interp1(x', I, time);

for i = 1 : length(I2)

    if isnan(I2(i))
        I2(i) = I2(i-1);
    end
    
end


%% entities declaration
A  = (-(1/M)*mpump)-1/(M*c*R);
B  = (1/M)*mpump;
C  = 1;
Bd = [1/(M*c*R) Sf/M*c Sf/M*c Sf/M*c Sf/M*c Q/M*c];
B = [B Bd];
D = 0;
model = ss(A,B,C,D);

%% controller definition

model = setmpcsignals(model, 'MV', 1, 'MD', 2:7);
mpc1v2 = mpc(model, Ts, P, m);
% load('mpc1v2.mat');
mpc1v2.Weights.OutputVariables = 1;
mpc1v2.Weights.ECR = 1000;
mpc1v2.Weights.ManipulatedVariables = 0.05;
mpc1v2.Weights.ManipulatedVariablesRate = 0.1;

mpc1v2.ManipulatedVariables.Min = 0;
mpc1v2.ManipulatedVariables.RateMin = -4;
mpc1v2.ManipulatedVariables.Max = 50; % A voir si on peut tolrere une telle température
mpc1v2.ManipulatedVariables.RateMax = 4;
mpc1v2.ManipulatedVariables.ScaleFactor = 60;

mpc1v2.OutputVariables.Min = 10;
mpc1v2.OutputVariables.MinECR = 1;
mpc1v2.OutputVariables.Max = 23.5;
mpc1v2.OutputVariables.MaxECR = 2;
mpc1v2.OutputVariables.ScaleFactor = 80;

mpc1v2.DisturbanceVariables(1).ScaleFactor = 80;
mpc1v2.DisturbanceVariables(2).ScaleFactor = 200;
mpc1v2.DisturbanceVariables(3).ScaleFactor = 200;
mpc1v2.DisturbanceVariables(4).ScaleFactor = 200;
mpc1v2.DisturbanceVariables(5).ScaleFactor = 200;
mpc1v2.DisturbanceVariables(5).ScaleFactor = 300;


params = mpcsimopt(mpc1v2);
params.MDLookAhead = 'on';
params.RefLookAhead = 'on';

v1 = [zeros(length(time), 6)];
v2 = [AmbianTemperatureInterpol' I2' I2' I2' I2' N'];
[y1, t1, u1] = sim(mpc1v2, Tstop/Ts, r',v1, params);
[y2, t2, u2] = sim(mpc1v2, Tstop/Ts, r',v2, params);

%% plot résults
figure(1);
plot(t1, r, 'green')
grid minor 
hold on
plot(t1, y1, 'red')
hold on 
stairs(t1, u1, 'blue')
hold on 
plot(t1, AmbianTemperatureInterpol, 'black')
title('temperature zone, ambiante, chauffage sans disturbances')


figure(2);
plot(t2, r, 'green')
grid minor 
hold on
plot(t2, y2, 'red')
hold on 
stairs(t2, u2, 'blue')
hold on 
plot(t2, AmbianTemperatureInterpol, 'black')
title('temperature zone, ambiante, chauffage avec disturbances')

%% Calcul de la consommation énergétique par le chauffage
Qheat1 = mpump*cpf*u1 - AmbianTemperatureInterpol';
Qheat2 = mpump*cpf*u2 - AmbianTemperatureInterpol';

figure;
plot(time, Qheat1/1000, '-.', 'LineWidth', 1.5) % on divise Qheat1 par 1000 pour la ploter en kw
hold on
plot(time, Qheat2/1000, 'LineWidth', 1.1) % on divise Qheat2 par 1000 pour la ploter en kw
grid minor
legend('Puissance sans disturbances Radiation solaire, temperature exterieure et activité humaine', 'avec activité humaine, météo');
title('Puissance instantannée consommée par le chauffage')

% La puissance consommée par jour pour le chauffage est à peu pres 18 kw
% contre 24 kw pour Trnsys

% 
% %% PID Trnsys for comparison