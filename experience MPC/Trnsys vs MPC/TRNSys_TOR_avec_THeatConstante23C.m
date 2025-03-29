clear all
clc
% close all

%% load file with gain

ville = 'Londre';
filePath = 'ResultsWithGainTheatOn.txt';
f = readtable(filePath);
PvProduction = f(:, 2); PvProduction = table2array(PvProduction);
HeatingPower = f(:, 3); HeatingPower = table2array(HeatingPower);
SolaRadiation = f(:, 4); SolaRadiation = table2array(SolaRadiation);
AmbianTemperature = f(:, 5); AmbianTemperature = table2array(AmbianTemperature);
ZoneTemperature = f(:, 6); ZoneTemperature = table2array(ZoneTemperature);
HeaterTemperature = f(:, 7); HeaterTemperature = table2array(HeaterTemperature);

%% load file without gain
filePath2 = 'ResultsWithoutGainTheatOn.txt';
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

timePlot = (0:119);

% With gain
AmbianTemperature = AmbianTemperature(1:120);
ZoneTemperature = ZoneTemperature(1:120);
HeaterTemperature = HeaterTemperature(1:120);
SolaRadiation = SolaRadiation(1:120);
PvProduction = PvProduction(1:120);
HeatingPower = HeatingPower(1:120);

% Without gain
AmbianTemperature2 = AmbianTemperature2(1:120);
ZoneTemperature2 = ZoneTemperature2(1:120);
HeaterTemperature2 = HeaterTemperature2(1:120);
SolaRadiation2 = SolaRadiation2(1:120);
PvProduction2 = PvProduction2(1:120);
HeatingPower2 = HeatingPower2(1:120);

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

%% référence (à refaire)
r = [zeros(1, 9), ones(1, 4), zeros(1, 1), ones(1, 5), zeros(1, 5)];
r = repmat(r, 1, 5);

r = 23 *r;

%% Plot figure 1
figure
subplot(2, 1, 1)
plot(timePlot, AmbianTemperature, 'color', 'black');
hold on
plot(timePlot, ZoneTemperature, 'color', 'blue')
hold on 
plot(timePlot, HeaterTemperature, 'color', 'red')
hold on 
plot(timePlot, r, 'color', 'green')
legend('reférence', 'Temperatture ambiante', 'Temperature zone', 'Temperature heater', 'FontSize', 14)
xlabel({'Time', '(h)'})
ylabel({'Temperature','C°'})
set(gca,'FontSize',14)
title(['Temperature ambiante,  de la zone et du heater en mode ON controlé par un simple thermostat pour la ville de ', ville]);
grid on

subplot(2, 1, 2)
plot(timePlot, SolaRadiation, 'color', 'blue', 'LineWidth', 2);
hold on 
plot(timePlot, PvProduction, 'color', 'green', 'LineWidth', 1.3)
hold on
plot(timePlot, HeatingPower, '-.r')
grid on
legend('radiation solaire', 'Production pv', 'puissance consommée par le heater', 'FontSize',14)
xlabel({'Time', '(h)'})
ylabel({'Puissance', '(kw)'})
set(gca,'FontSize',14)
title(['radiation solaire, production pv et consommation énergétique annuelle pour la ville de ', ville, ' avec une activité humaine']);

%% figure 2
figure
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

% %% calcul de la consommation energetique a partir de la puissance
% E1 = cumsum(HeatingPower);
% E2 = cumsum(HeatingPower2);
% figure(3)
% plot(timePlot, E1, 'blue')
% hold on 
% plot(timePlot, E2, 'red')
% hold on
% plot(timePlot, E1(end, 1) * linspace(1, 1, length(E1)), 'color', 'blue','LineStyle', '--')
% hold on
% plot(timePlot,E2(end, 1) * linspace(1, 1, length(E2)), 'color', 'red','LineStyle', '--')
% legend('Heater Energy consumption with humains gain', 'Heater Energy consumption without humans activity')
% % xticks([0 720 1441 2162 2881 3600 4321 5040 5760 6481 7200 7922 length(time)])
% % xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% grid minor
% title('difference entre consommation du mois d''aout with gain and without gain')

