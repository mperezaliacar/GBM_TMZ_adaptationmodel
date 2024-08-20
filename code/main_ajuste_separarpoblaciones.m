clear all
close all
clc

rng(0); % seed (reproducibility)

global xref xref_repair N days

%%%%%%%%%%%%%%%%%%% Exp data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
cd ..
load expDATA/expDATA_ind
name = 'paperTea';
cd ajuste/Nivel1

N = size(datos_tratados,1);
days = time;

filename1 = 'results_modelselection/results_lsqnonlin_f1_paperTea5.mat';

%% FASE 2

xref(1) = 142; %L
xref(2) = 2.63; %n
xref(3) = 1e-3; %kappa
xref(4) = 0.34; %dth
xref(5) = 0.45; %DeltaD
xref(6) = 5.8e9; %tauR
xref(7) = 1.04; %H

xref_repair(1) = 6.9e-4; %lambda
xref_repair(2) = 1; %beta

%%
LB = zeros(1,length(xref));
UB = inf*ones(1,length(xref));
X0 = ones(1,length(xref)); 

stol = 1e-10;
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','StepTolerance',stol);

FUN = @(x) objfun_f2_ind(x,filename1,datos_tratados);

[x,resnorm,fval] = lsqnonlin(FUN,X0,LB,UB,options);

parameters = x.*xref;


filename2 = strcat('results_modelselection/parameters_f2_ind');
save(filename2,'parameters')    




