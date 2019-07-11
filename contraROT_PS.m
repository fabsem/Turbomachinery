clear all
close all
clc

%% Add Paths
%addPaths;
addpath('liebleinCorrelations/')
addpath('howellCorrelations/')
addpath('losses/')
addpath('IGV/')
addpath('traupel/')
addpath('blades/')

%% Optimization
%work rotRatio n Mw1_tip Re1 Re2
lb = [0.3 0.8 2600 0.7 6e5 6e5];
ub = [0.7 1.5 3800 0.8 1.5e6 1.5e6];

fun = @contraROTTO;
nvars = 6;

%Starts ParallelPool for Parallel Computation
%parpool;

%Max simulation time
days = 7;
MaxTime = 60 * 60 * 24 * days;

options = optimoptions('particleswarm','SwarmSize',100,'Display','Iter',...
                       'FunctionTolerance',1e-3,'MaxIterations',200*nvars,...
                       'MaxTime',MaxTime,'PlotFcn','pswplotbestf',...
                       'UseParallel',true,'UseVectorized',false);
                                      
[optimalSolution,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);