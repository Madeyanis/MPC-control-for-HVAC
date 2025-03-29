clear all
clc
close all

%% load file with gain

ville = 'Londre';
filePath = 'ResultsWithGain.txt';
f = readtable(filePath);
PvProduction = f(:, 2); PvProduction = table2array(PvProduction);
HeatingPower = f(:, 3); HeatingPower = table2array(HeatingPower);
SolaRadiation = f(:, 4); SolaRadiation = table2array(SolaRadiation);
AmbianTemperature = f(:, 5); AmbianTemperature = table2array(AmbianTemperature);
ZoneTemperature = f(:, 6); ZoneTemperature = table2array(ZoneTemperature);
HeaterTemperature = f(:, 7); HeaterTemperature = table2array(HeaterTemperature);

%% load file without gain
filePath2 = 'ResultsWithoutGain.txt';
f2 = readtable(filePath2);
PvProduction2 = f2(:, 2); PvProduction2 = table2array(PvProduction2);
HeatingPower2 = f2(:, 3); HeatingPower2 = table2array(HeatingPower2);
SolaRadiation2 = f2(:, 4); SolaRadiation2 = table2array(SolaRadiation2);
AmbianTemperature2 = f2(:, 5); AmbianTemperature2 = table2array(AmbianTemperature2);
ZoneTemperature2 = f2(:, 6); ZoneTemperature2 = table2array(ZoneTemperature2);
HeaterTemperature2 = f2(:, 7); HeaterTemperature2 = table2array(HeaterTemperature2);

%% Some manipulations for changing time to predict
keySet = {'day', 'month', 'year'};
valueSet = [30, 720, 8760];
TimeToPredict = containers.Map(keySet,valueSet);
choix = 'year';


%% Choix de la pèriode de l'année personnalisé (remarque qu'il faut appliquer soit ça ou les periodes avec dictionnaire)
% on part sur la partie ou le gain fait la différence, (faut changer
% le timetoplot)
% prendre le moi de septembre ce qui correspond à 6552:744

timePlot = (5112:(5112+744));

% With gain
AmbianTemperature = AmbianTemperature(5112:(5112+744));
ZoneTemperature = ZoneTemperature(5112:(5112+744));
HeaterTemperature = HeaterTemperature(5112:(5112+744));
SolaRadiation = SolaRadiation(5112:(5112+744));
PvProduction = PvProduction(5112:(5112+744));
HeatingPower = HeatingPower(5112:(5112+744));

% Without gain
AmbianTemperature2 = AmbianTemperature2(5112:(5112+744));
ZoneTemperature2 = ZoneTemperature2(5112:(5112+744));
HeaterTemperature2 = HeaterTemperature2(5112:(5112+744));
SolaRadiation2 = SolaRadiation2(5112:(5112+744));
PvProduction2 = PvProduction2(5112:(5112+744));
HeatingPower2 = HeatingPower2(5112:(5112+744));

% Periode de 1 jours ou 1 mois ou 1 année
% timePlot = (0 : (TimeToPredict(choix)-1));
% % With gain
% AmbianTemperature = AmbianTemperature(1:TimeToPredict(choix));
% ZoneTemperature = ZoneTemperature(1:TimeToPredict(choix));
% HeaterTemperature = HeaterTemperature(1:TimeToPredict(choix));
% 
% SolaRadiation = SolaRadiation(1:TimeToPredict(choix));
% PvProduction = PvProduction(1:TimeToPredict(choix));
% HeatingPower = HeatingPower(1:TimeToPredict(choix));
% 
% % Without gain
% AmbianTemperature2 = AmbianTemperature2(1:TimeToPredict(choix));
% ZoneTemperature2 = ZoneTemperature2(1:TimeToPredict(choix));
% HeaterTemperature2 = HeaterTemperature2(1:TimeToPredict(choix));
% 
% SolaRadiation2 = SolaRadiation2(1:TimeToPredict(choix));
% PvProduction2 = PvProduction2(1:TimeToPredict(choix));
% HeatingPower2 = HeatingPower2(1:TimeToPredict(choix));

%% Plot figure 1
figure(1)
subplot(2, 1, 1)
plot(timePlot, AmbianTemperature, 'color', 'black');
hold on
plot(timePlot, ZoneTemperature, 'color', 'blue')
hold on 
plot(timePlot, HeaterTemperature, 'color', 'red')
legend('temperatture ambiante', 'Temperature zone', 'Temperature heater')
title(['Temperature ambiante,  de la zone et du heater en mode OFF controlé par un simple thermostat pour la ville de ', ville]);
grid on
subplot(2, 1, 2)
plot(timePlot, SolaRadiation, 'color', 'blue', 'LineWidth', 2);
hold on 
plot(timePlot, PvProduction, 'color', 'green', 'LineWidth', 1.3)
hold on
plot(timePlot, HeatingPower, '-.r')
grid on
legend('radiation solaire', 'Production pv', 'puissance consommée par le heater')
title(['radiation solaire, production pv et consommation énergétique annuelle pour la ville de ', ville, ' avec une activité humaine']);

%% figure 2
figure(2)
subplot(2, 1, 1)
plot(timePlot, AmbianTemperature2, 'color', 'black');
hold on
plot(timePlot, ZoneTemperature2, 'color', 'blue')
hold on 
plot(timePlot, HeaterTemperature2, 'color', 'red')
legend('temperatture ambiante', 'Temperature zone', 'Temperature heater')
title(['Temperature ambiante,  de la zone et du heater en mode OFF controlé par un simple thermostat pour la ville de ', ville, ' sans activité humaine']);
grid on
subplot(2, 1, 2)
plot(timePlot, SolaRadiation2, 'color', 'blue', 'LineWidth', 2);
hold on 
plot(timePlot, PvProduction2, 'color', 'green', 'LineWidth', 1.3)
hold on
plot(timePlot, HeatingPower2, '-.r')
grid on
legend('radiation solaire', 'Production pv', 'puissance consommée par le heater')
title(['radiation solaire, production pv et consommation énergétique annuelle pour la ville de ', ville]);

%% calcul de la consommation energetique a partir de la puissance
E1 = cumsum(HeatingPower);
E2 = cumsum(HeatingPower2);
figure(3)
plot(E1, 'blue')
hold on 
plot(E2, 'red')
hold on
plot(E1(end, 1) * linspace(1, 1, length(E1)), 'color', 'blue','LineStyle', '--')
hold on
plot(E2(end, 1) * linspace(1, 1, length(E2)), 'color', 'red','LineStyle', '--')
legend('Heater Energy consumption with humains gain', 'Heater Energy consumption without humans activity')
% xticks([0 720 1441 2162 2881 3600 4321 5040 5760 6481 7200 7922 length(time)])
% xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
grid minor
title('difference entre consommation du mois d''aout with gain and without gain')


%% Building definition
Volume = 84;
DensiteAir = 1.204;
% M = Masse d'air à l'interieur de la zone (kg)
M = Volume * DensiteAir;
% c = Capacité calorifique de l'air à pression constante (j/kg.K)
c = 1256;
% m = Débit du fluid entrant dans le heater
mpump = 50;
% Theat = température du heater (commande à réguler)
% Q = Apport d'un seul occupant
Q = 150;
% N = Nombre d'occupant
N = 2;
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
% Chauffage
mcpf = 4.19 * 1000;

% %% Some variables for Predictions
% Ts = 0.5;
% P = 24;
% m = 8;
% % Tstop permet de définir la longueur de temps qu'on veut prédire
% Tstop = TimeToPredict(choix);
% 
% time = (0 : Ts : (Tstop)-Ts);
% 
% % Interpollser par rapport à la longeur souhaitée
% x = 0:(TimeToPredict(choix)-1); % x = longeur souhaitée - 1
% AmbianTemperature = interp1(x', AmbianTemperature, time);
% 
% for i = 1 : length(AmbianTemperature)
% 
%     if isnan(AmbianTemperature(i))
%         AmbianTemperature(i) = AmbianTemperature(i-1);
%     end
%     
% end
% 
% SolaRadiation = interp1(x', SolaRadiation, time);
% 
% for i = 1 : length(SolaRadiation)
% 
%     if isnan(SolaRadiation(i))
%         SolaRadiation(i) = SolaRadiation(i-1);
%     end
%     
% end
% 
% %% Reference definition
% r = 21*[zeros(1, 15) ones(1, 9) zeros(1, 1) ones(1, 11) zeros(1, 12)];
% % plot(r)
% rr = repmat(r, 1, 30);
% % plot(rr)

%% Problem definition 
% dx/dt = Ax + Bu + Bd
% y = Cx
% x = [Tz] state représente la température de la zone
% y = [Tz] sortie     //     //    //       // //  //
% u = [Theat] control signal: température du chauffage
% d = [TOut Ie Iw Is In N] les entrées non commandables mais mesurables


%% entities declaration
% A  = (-(1/M)*mpump)-1/(M*c*R);
% B  = (1/M)*mpump;
% C  = 1;
% Bd = [1/(M*c*R) Sf/M*c Sf/M*c Sf/M*c Sf/M*c Q/M*c];
% B = [B Bd];
% D = 0;
% model = ss(A,B,C,D);


%% Partie prédiction par MPC
% 
% model = setmpcsignals(model, 'MV', 1, 'MD', 2:7);
% mpc1v2 = mpc(model, Ts, 4, 2);
% % load('mpc1v2.mat');
% mpc1v2.Weights.OutputVariables = 1;
% mpc1v2.Weights.ECR = 1000;
% mpc1v2.Weights.ManipulatedVariables = 0.005;
% mpc1v2.Weights.ManipulatedVariablesRate = 0.1;
% 
% mpc1v2.ManipulatedVariables.Min = 0;
% mpc1v2.ManipulatedVariables.RateMin = -10;
% mpc1v2.ManipulatedVariables.Max = 50; % A voir si on peut tolrere une telle température
% mpc1v2.ManipulatedVariables.RateMax = 10;
% mpc1v2.ManipulatedVariables.ScaleFactor = 60;
% 
% mpc1v2.OutputVariables.Min = 10;
% mpc1v2.OutputVariables.MinECR = 1;
% mpc1v2.OutputVariables.Max = 21.5;
% mpc1v2.OutputVariables.MaxECR = 2;
% mpc1v2.OutputVariables.ScaleFactor = 80;
% 
% mpc1v2.DisturbanceVariables(1).ScaleFactor = 80;
% mpc1v2.DisturbanceVariables(2).ScaleFactor = 700;
% mpc1v2.DisturbanceVariables(3).ScaleFactor = 700;
% mpc1v2.DisturbanceVariables(4).ScaleFactor = 700;
% mpc1v2.DisturbanceVariables(5).ScaleFactor = 700;
% mpc1v2.DisturbanceVariables(5).ScaleFactor = 250;
% 
% 
% params = mpcsimopt(mpc1v2);
% params.MDLookAhead = 'on';
% params.RefLookAhead = 'on';
% 
% v = [AmbianTemperature' zeros(length(time), 5)];
% 
% [y, t, u] = sim(mpc1v2, Tstop/Ts, rr',v, params);
% % 
% % %% Consommation energétique
% Qheat = ((mcpf*mpump) * u) - AmbianTemperature;
% Qheat = Qheat * 0.00000027777777778;
% %% plot résults
% 
% plot(t, rr, 'green')
% grid minor 
% hold on
% plot(t, y, 'red')
% hold on 
% stairs(t, u, 'blue')
% hold on 
% plot(t, AmbianTemperature, 'black')
% figure
% plot(Qheat, 'magenta')
% title ('Consommation energétique')
% grid minor
% 
% % 
% % %% PID Trnsys for comparison