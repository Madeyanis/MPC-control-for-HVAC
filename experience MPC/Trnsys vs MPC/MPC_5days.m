clear all
clc
close all


%% Cette experience est realisee pour une duree de 5 jours


%% Problem definition 
% dx/dt = Ax + Bu + Bd
% y = Cx
% x = [Tz] state représente la température de la zone
% y = [Tz] sortie     //     //    //       // //  //
% u = [Theat] control signal: température du chauffage
% d = [TOut Ie Iw Is In N] les entrées non commandables mais mesurables

%% Load fichier
data = readtable('ResultsWithoutGainTheatOn.txt');
data = table2array(data);
data(1, :) = [];

%% Nombre de journée a prévoir
jour = 5;

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
N = repmat(N, 1, jour);
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
R = R*0.8;
% Ppv : Production panneaux photovoltaiques 
PvProduction = data(:, 2); PvProduction = PvProduction * 1000;
% Radiation solaire 
SolaRadiation = data(:, 4); SolaRadiation = SolaRadiation * 1000;
% temperature ambiante
AmbianTemperature = data(:, 5);
% Chauffage
cpf = (4.19/3.6);

%% Some variables for Mpc
Ts = 0.2;
P = 5;
m = 3;
Tstop = 24*jour;

time = (0 : Ts : (Tstop)-Ts);

% Reference definition
r = N/2;
r = 19 *r;
% r(length(time)+1:end)= [];

%% time manipulation for variables

AmbianTemperature = AmbianTemperature(1:(24*jour));
PvProduction = PvProduction(1:(24*jour));
SolaRadiation = SolaRadiation(1:(24*jour));




x = 0:((24*jour)-1);
AmbianTemperature = interp1(x', AmbianTemperature, time);
SolaRadiation = interp1(x', SolaRadiation, time);
PvProduction = interp1(x', PvProduction, time);

for i = 1 : length(AmbianTemperature)

    if isnan(AmbianTemperature(i))
        AmbianTemperature(i) = AmbianTemperature(i-1);
    end
    if isnan(SolaRadiation(i))
        SolaRadiation(i) = SolaRadiation(i-1);
    end
    if isnan(PvProduction(i))
        PvProduction(i) = PvProduction(i-1);
    end
    
end


%% Normaliser les unités
PvProduction = PvProduction / 1000; %% Pour avoir les KW

%% entities declaration
A  = (-(1/M)*mpump)-1/(M*c*R);
B  = (1/M)*mpump;
C  = 1;
Bd = [1/(M*c*R) Sf/M*c Sf/M*c Sf/M*c Sf/M*c Q/M*c];
B = [B Bd];
D = zeros(1, 7);
model = ss(A,B,C,D);

%% controller definition

temperatureMaxDuHeater = 34;

model = setmpcsignals(model, 'MV', 1, 'MD', 2:7);
mpc1v2 = mpc(model, Ts, P, m);
% load('mpc1v2.mat');
mpc1v2.Weights.OutputVariables = 10;
mpc1v2.Weights.ECR = 1000;
mpc1v2.Weights.ManipulatedVariables = 0.1;
mpc1v2.Weights.ManipulatedVariablesRate = 0.2;

mpc1v2.ManipulatedVariables.Min = AmbianTemperature;
% mpc1v2.ManipulatedVariables.RateMin = -(temperatureMaxDuHeater/5);
mpc1v2.ManipulatedVariables.RateMin = -6;
mpc1v2.ManipulatedVariables.Max = temperatureMaxDuHeater; % A voir si on peut tolrere une telle température
% mpc1v2.ManipulatedVariables.RateMax = (temperatureMaxDuHeater/5);
mpc1v2.ManipulatedVariables.RateMax = 6;
mpc1v2.ManipulatedVariables.ScaleFactor = 60;

mpc1v2.OutputVariables.Min = -Inf;
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
v2 = [AmbianTemperature' SolaRadiation' SolaRadiation' SolaRadiation' SolaRadiation' N'];
[Tzone1, t1, Theat1] = sim(mpc1v2, Tstop/Ts, r',v1, params);
[Tzone2, t2, Theat2] = sim(mpc1v2, Tstop/Ts, r',v2, params);

%% plot résults
figure;
subplot(211)
plot(t1, r, 'green')
grid minor 
hold on
plot(t1, Tzone1, 'red')
hold on 
stairs(t1, Theat1, 'blue')
hold on 
plot(t1, AmbianTemperature, 'black')
legend('Référence', 'Temperature de la zone', 'Temperature du heater', 'Temperature ambiante', 'FontSize', 12)
xlabel({'Time', '(h)'})
ylabel({'Temperature','C°'})
set(gca,'FontSize',14)
ylim([0 (temperatureMaxDuHeater + 5)])
title('temperature zone, ambiante, chauffage sans disturbances')

subplot(212)
plot(t2, r, 'green')
grid minor 
hold on
plot(t2, Tzone2, 'red')
hold on 
stairs(t2, Theat2, 'blue')
hold on 
plot(t2, AmbianTemperature, 'black')
legend('Référence', 'Temperature de la zone', 'Temperature du heater', 'Temperature ambiante', 'FontSize', 12)
xlabel({'Time', '(h)'})
ylabel({'Temperature','C°'})
set(gca,'FontSize',14)
ylim([0 (temperatureMaxDuHeater + 5)])
title('temperature zone, ambiante, chauffage avec disturbances')

%% Calcul de la consommation énergétique par le chauffage
Qheat1 = mpump*cpf*(Theat1 - AmbianTemperature');
Qheat2 = mpump*cpf*(Theat2 - AmbianTemperature');

for i = 1 : length(Qheat1)
    if Qheat1(i) < 0
        Qheat1(i) = 0;
    end
    if Qheat2(i) < 0 
        Qheat2(i) = 0;
    end
end
        

figure;
plot(time, Qheat1/1000, '-.', 'LineWidth', 1.1) % on divise Qheat1 par 1000 pour la ploter en kw
hold on
plot(time, Qheat2/1000, '-.', 'LineWidth', 1.1) % on divise Qheat2 par 1000 pour la ploter en kw
hold on 
plot(time, PvProduction, 'green', 'LineWidth', 1.5)
grid minor
legend('Power consumed without disturbances', 'Power consumed with disturbances', 'Pv Prodcution', 'FontSize', 14);
xlabel({'Time', '(h)'})
ylabel({'Puissance','(Kw)'})
set(gca,'FontSize',14)
title('Puissance instantannée consommée par le chauffage')



Eheat1 = sum(Qheat1(15:end))/(5*1000);
Eheat2 = sum(Qheat2(15:end))/(5*1000);


%% plot résults
figure;
subplot(211)
plot(t2, r, 'green')
grid minor 
hold on
plot(t2, Tzone2, 'red')
hold on 
stairs(t2, Theat2, 'blue')
hold on 
plot(t2, AmbianTemperature, 'black')
legend('Référence', 'Temperature de la zone', 'Temperature du heater', 'Temperature ambiante', 'FontSize', 12)
xlabel({'Time', '(h)'})
ylabel({'Temperature','C°'})
set(gca,'FontSize',14)
ylim([0 (temperatureMaxDuHeater + 5)])
title('temperature zone, ambiante, chauffage avec disturbances')

subplot(212)
plot(time, Qheat2/1000, '-.', 'LineWidth', 1.1) % on divise Qheat2 par 1000 pour la ploter en kw
hold on 
plot(time, PvProduction, 'green', 'LineWidth', 1.5)
grid minor
legend('Heater power consumed', 'Pv Prodcution', 'FontSize', 14);
xlabel({'Time', '(h)'})
ylabel({'Puissance','(Kw)'})
set(gca,'FontSize',14)
title('Puissance instantannée consommée par le chauffage et production PV')

%% definition of second criteria 
% Comfort criteria

% j1p1(1, :) = r(40:(11.8*5));
% j1p1(2, :) = Tzone2(8*5:11.8*5);
% j1p2(1, :) = r(13.2*5:18*5);
% j1p2(2, :) = Tzone2(13.2*5:18*5);
% 
% figure;
% plot(j1p2(1, :))
% hold on 
% plot(j1p2(2, :))
% e1 = abs(sum(j1p1(1, :)-j1p1(2, :)))
% e2 = abs(sum(j1p2(1, :)-j1p2(2, :)))