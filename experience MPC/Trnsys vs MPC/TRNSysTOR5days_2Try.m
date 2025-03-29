clear all
clc
close all


%% load file with gain

ville = 'Londre';
filePath = 'TOR3.txt';
f = readtable(filePath);

AmbianTemperature = f(:, 5); AmbianTemperature = table2array(AmbianTemperature);
ZoneTemperature = f(:, 6); ZoneTemperature = table2array(ZoneTemperature);
HeaterTemperature = f(:, 7); HeaterTemperature = table2array(HeaterTemperature);
timeToPlot = f(:, 1); timeToPlot = table2array(timeToPlot);
PvProduction = f(:, 2); PvProduction = table2array(PvProduction);
HeatingPower = f(:, 3); HeatingPower = table2array(HeatingPower);


%% référence (à refaire parce qu'il faut faire la réference sur 2 fois plus de samples)
r = [0 zeros(1, 39), ones(1, 21), zeros(1, 5), ones(1, 25), zeros(1, 29)];
r = repmat(r, 1, 5);
r = [r zeros(1, 1)];
r = 19 *r';

%%
% 
% x = 0:3600;
% ZoneTemperature = interp1(timeToPlot, ZoneTemperature, x);

%% Plot figure 1
figure;
subplot(211)
plot(timeToPlot, AmbianTemperature, 'black')
hold on 
plot(timeToPlot, HeaterTemperature, 'red')
hold on
plot(timeToPlot, ZoneTemperature, 'blue')
hold on
plot(timeToPlot, r, 'green')
grid minor
legend('Tout', 'Theat', 'Tzone')
xlabel('t(h)')
ylabel({'Temperature','C°'})
ylim([0 35])
set(gca,'FontSize',14)
title('Temperatures with TOR Control')

subplot(212)
plot(timeToPlot, HeatingPower, 'red')
hold on 
plot(timeToPlot, PvProduction, 'green', 'LineWidth', 1.6)
legend('Heater Consumption','Production PV')
grid minor
xlabel('t(h)')
ylabel({'Power','Kw'})
ylim([0 2])
set(gca,'FontSize',14)
title('heater energy with TOR Control and pv prodution')

E = 0;
E = sum(HeatingPower/5);

%% Comfort criteria

j1p1(1, :) = r(8.2*5:12.2*5);
j1p1(2, :) = ZoneTemperature(8.2*5:12.2*5);
j1p2(1, :) = r(13.4*5:18.2*5);
j1p2(2, :) = ZoneTemperature(13.4*5:18.2*5);


j2p1(1, :) = r(32.2*5:36.2*5);
j2p1(2, :) = ZoneTemperature(32.2*5:36.2*5);
j2p2(1, :) = r(37.4*5:42.2*5);
j2p2(2, :) = ZoneTemperature(37.4*5:42.2*5);

figure;
plot(j1p2(1, :))
hold on 
plot(j1p2(2, :))
e1 = abs(sum(j1p1(1, :)-j1p1(2, :)));
e2 = abs(sum(j1p2(1, :)-j1p2(2, :)));

figure;
plot(j2p2(1, :))
hold on 
plot(j2p2(2, :))
e3 = abs(sum(j2p1(1, :)-j2p1(2, :)));
e4 = abs(sum(j2p2(1, :)-j2p2(2, :)));
e = e1 + e2 + e3 + e4