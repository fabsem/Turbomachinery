clear all
close all
clc

%% Add Paths
addpath('../liebleinCorrelations/')
addpath('../howellCorrelations/')
addpath('../losses/')
addpath('../IGV/')
addpath('../traupel/')
addpath('../blades/')
addpath('../velocityTriangles')
addpath('../PSO')
addpath('../dataAssignments/')

%% Optimization
%work rotRatio n Mw1_tip Re1 Re2 alfa1
lb = [0.3 0.7 2600 0.7 6e5 6e5 -5  -5  -5];
ub = [0.7 1.5 3800 0.88 2e6 2e6  7  7  7];

fun = @contraROT_objFUN;
nvars = 9;

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
