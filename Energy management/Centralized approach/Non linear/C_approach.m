clear all
close all
clc


nx = 3;
nu = 6;
ny = 3;
nlobj = nlmpc(nx, ny, nu);

nlobj.Model.StateFcn = @StateJcbFcnBc;

dt = 0.1;
p = 24;
nlobj.Ts = dt;
nlobj.PredictionHorizon = p;
nlobj.ControlHorizon = round(p/8);
nlobj.Optimization.ReplaceStandardCost = true;

nlobj.Model.OutputFcn = @myOutputFunction;



nlobj.Optimization.CustomCostFcn = @(X, U, data) 10*sum(sum((U(1:end-1, data.MVIndex(1)).^2)));
