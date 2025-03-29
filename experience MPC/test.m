clear all
close all
clc


%% With Humains Gain
filePath = 'ResultsWithGain.txt';
f = readtable(filePath);
PvEnergy = f(:, 2); PvEnergy = table2array(PvEnergy);
HeatingEnergy = f(:, 3); HeatingEnergy = table2array(HeatingEnergy);
SolaRadiation = f(:, 4); SolaRadiation = table2array(SolaRadiation);
AmbianTemperature = f(:, 5); AmbianTemperature = table2array(AmbianTemperature);
ZoneTemperature = f(:, 6); ZoneTemperature = table2array(ZoneTemperature);
HeaterTemperature = f(:, 7); HeaterTemperature = table2array(HeaterTemperature);
disp(['heating energy with gain: ', num2str(HeatingEnergy(end, 1)), ' kwh'])

%% Without humans Activity
filePath = 'ResultsWithoutGain.txt';
f = readtable(filePath);
PvEnergy2 = f(:, 2); PvEnergy2 = table2array(PvEnergy2);
HeatingEnergy2 = f(:, 3); HeatingEnergy2 = table2array(HeatingEnergy2);
SolaRadiation2 = f(:, 4); SolaRadiation2 = table2array(SolaRadiation2);
AmbianTemperature2 = f(:, 5); AmbianTemperature2 = table2array(AmbianTemperature2);
ZoneTemperature2 = f(:, 6); ZoneTemperature2 = table2array(ZoneTemperature2);
HeaterTemperature2 = f(:, 7); HeaterTemperature2 = table2array(HeaterTemperature2);
disp(['heating energy without activity: ', num2str(HeatingEnergy2(end, 1)), ' kwh'])

%% Energy comparaison
figure(1)
plot(HeatingEnergy)
hold on 
plot(HeatingEnergy2)
hold on
plot(HeatingEnergy(end, 1) * linspace(1, 1, length(HeatingEnergy)), 'color', 'red','LineStyle', '--')
hold on
plot(HeatingEnergy2(end, 1) * linspace(1, 1, length(HeatingEnergy2)), 'color', 'red','LineStyle', '--')
grid minor
legend('Energy with humains gain', 'Energy without humans activity')
xticks([0 720 1441 2162 2881 3600 4321 5040 5760 6481 7200 7922 length(PvEnergy)])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
title('Heater consumption with and without humains activity')

% %% Zone Temperature
% figure(2)
% plot(ZoneTemperature)
% hold on 
% plot(ZoneTemperature2)
% grid minor
% legend('Zone temperature with humains gain', 'Zone Temperature without humans activity')
% xticks([0 720 1441 2162 2881 3600 4321 5040 5760 6481 7200 7922 length(PvEnergy)])
% xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12'})
% title('Zone temperature with and without humains activity without Heater')

