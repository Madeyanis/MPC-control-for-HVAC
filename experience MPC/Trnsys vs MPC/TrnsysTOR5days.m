clear all
clc
close all


%% load file with gain

ville = 'Londre';
filePath = 'TOR5daysWithGain.txt';
f = readtable(filePath);

AmbianTemperature = f(:, 5); AmbianTemperature = table2array(AmbianTemperature);
ZoneTemperature = f(:, 4); ZoneTemperature = table2array(ZoneTemperature);
HeaterTemperature = f(:, 2); HeaterTemperature = table2array(HeaterTemperature);
timeToPlot = f(:, 1); timeToPlot = table2array(timeToPlot);
PvProduction = f(:, 6); PvProduction = table2array(PvProduction);
HeatingPower = f(:, 7); HeatingPower = table2array(HeatingPower);


%% référence (à refaire parce qu'il faut faire la réference sur 2 fois plus de samples)
r = [zeros(1, 8*50), ones(1, 4*50), zeros(1, 1*50), ones(1, 5*50), zeros(1, (6*50))];
r = repmat(r, 1, 5);
r = [r zeros(1, 1)];
r = 19 *r;

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
for i = 1 :50: length(HeatingPower)
    E = E + HeatingPower(i);
end

mmm = zeros(1, 120);

for i = 1 : 119
    mmm(i) = mean(HeaterTemperature((i*50+1):(i*50+50)));
end

    