clear all
close all
clc

ts = 1;

clear all
close all
clc

%% Declaration et construction des signaux PV et Pl
ts = 1;
a = 100;
b = 600;
pl = a +(b-a)*rand(1, 48);
pv = [zeros(1, 6) 250*ones(1, 4) 600*ones(1, 4) 250*ones(1, 4) zeros(1, 6)]; 
pv = [pv pv];


%% Construction state model
A = 1;
B = [0 ts];
C = [1; 0];
D = [0 0; 1 1];
model = ss(A, B, C, D);
model.InputName = {'Pg', 'Pb'};
model.OutputName = {'Soc', 'y2'};
model.InputDelay = 0.001;


r1 = 500*ones(1,48)';
r2 = (pl - pv)';
r = [r1 r2];

%% create MPC controller object with sample time
mpc1_C = mpc(model, ts);
%% specify prediction horizon
mpc1_C.PredictionHorizon = 7;
%% specify control horizon
mpc1_C.ControlHorizon = 3;
%% specify nominal values for inputs and outputs
mpc1_C.Model.Nominal.U = [0;0];
mpc1_C.Model.Nominal.Y = [500;0];
%% specify scale factors for inputs and outputs
mpc1_C.MV(1).ScaleFactor = 1600;
mpc1_C.MV(2).ScaleFactor = 2000;
mpc1_C.OV(1).ScaleFactor = 600;
mpc1_C.OV(2).ScaleFactor = 1000;
%% specify constraints for MV and MV Rate
mpc1_C.MV(1).Min = -800;
mpc1_C.MV(1).Max = 800;
mpc1_C.MV(1).RateMin = -150;
mpc1_C.MV(1).RateMax = 150;
mpc1_C.MV(2).Min = -1000;
mpc1_C.MV(2).Max = 1000;
mpc1_C.MV(2).RateMin = -250;
mpc1_C.MV(2).RateMax = 250;
%% specify constraint softening for MV and MV Rate
mpc1_C.MV(1).MinECR = 1;
mpc1_C.MV(1).MaxECR = 1;
mpc1_C.MV(1).RateMinECR = 1;
mpc1_C.MV(1).RateMaxECR = 1;
mpc1_C.MV(2).MinECR = 1;
mpc1_C.MV(2).MaxECR = 1;
mpc1_C.MV(2).RateMinECR = 1;
mpc1_C.MV(2).RateMaxECR = 1;
%% specify constraints for OV
mpc1_C.OV(1).Min = 200;
mpc1_C.OV(1).Max = 800;
%% specify constraint softening for OV
mpc1_C.OV(1).MinECR = 0.1;
mpc1_C.OV(1).MaxECR = 0.1;
%% specify weights
mpc1_C.Weights.MV = [0 0];
mpc1_C.Weights.MVRate = [0.1 0.1];
mpc1_C.Weights.OV = [0.05 1];
mpc1_C.Weights.ECR = 100000;
%% specify simulation options
options = mpcsimopt();
options.RefLookAhead = 'on';
options.MDLookAhead = 'on';
options.Constraints = 'on';
options.OpenLoop = 'off';
%% run simulation
[y, t, u] = sim(mpc1_C, 48, r, options);

figure
subplot(211)
plot(t, y(:, 2))
grid minor
hold on
plot(t, r2)
hold on
plot(t, y(:, 1))
legend('Pb + pg', 'Pl - Pv', 'Soc')
xlabel('t(h)')
ylabel({'puissance','Kw'})
set(gca,'FontSize',14)
ylim([0 1000])
title('Outputs')

subplot(212)
plot(t, u(:, 1))
hold on
plot(t, u(:, 2))
legend('pg', 'pb')
ylim([-1000 1000])
grid minor
xlabel('t(h)')
ylabel({'puissance','Kw'})
set(gca,'FontSize',14)
title('Inputs')
