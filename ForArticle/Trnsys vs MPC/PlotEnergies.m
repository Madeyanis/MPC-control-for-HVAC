close all
clear all
clc


qheat = load('\\pers.stockage.univ-lorraine.fr\masdoua1\Bureau\Qheat1.mat');
heatP = load('\\pers.stockage.univ-lorraine.fr\masdoua1\Bureau\TOR_Power.mat');
qheat = cell2mat(struct2cell(qheat));
heatP = cell2mat(struct2cell(heatP));

t1 = 0:599;
t2 = 0:6000;

figure
subplot(211)
plot(t1, qheat/1000)
ylim([0 2])
grid minor
xlabel('t(h)')
ylabel({'Power','Kw'})
set(gca,'FontSize',14)
title('Heating power with MPC Controller')

subplot(212)
plot(t2, heatP)
ylim([0 2])
grid minor
xlabel('t(h)')
ylabel({'Power','Kw'})
set(gca,'FontSize',14)
title('Heating power with TOR Controller')


Eheat = 0;
for i = 1 : 5 :length(qheat)
    Eheat = Eheat + qheat(i);
end

E = 0;
for i = 1 :50: length(heatP)
    E = E + heatP(i);
end