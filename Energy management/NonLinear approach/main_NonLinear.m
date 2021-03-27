clear all
close all
clc


nx = 1;
nu = 2;
ny = 2;
nlobj = nlmpc(nx, ny, nu);

nlobj.Model.StateFcn = @StateFcnB;

dt = 0.1;
p = 24;
nlobj.Ts = dt;
nlobj.PredictionHorizon = p;
nlobj.ControlHorizon = round(p/8);
nlobj.Optimization.ReplaceStandardCost = true;

nlobj.Model.OutputFcn = @myOutputFunction;
x = 

for i = 1 : p
    
    y(i, :) = myOutputFunction(X(i, :)', U(i, :)');
    
end

nlobj.Optimization.CustomCostFcn = @(X, U, data) 10*sum(sum((U(1:end-1, data.MVIndex(1)).^2)));
