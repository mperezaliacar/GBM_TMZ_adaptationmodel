%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         optimisation_stage1.m                           %
%                                                                         %
%    This script performs the STAGE 1 of the optimisation of the model    %
%    parameters, that is, those related to spheroid growth, for           %
%    a given model candidate (defined by setting to zero the value of     %
%    some parameters) and computes the resulting values of the chi^2,     %
%    BIC and AIC that can be used to perform model selection.             %
%                                                                         %
%    The user must specify the model candidate as well as give a number   %
%    for each candidate (candidate_num) to identify the results.          %       
%                                                                         %
%    The results are stored in the results folder as:                     %
%    "results_optimisation_stage1_candidateX.mat" with X the number       %
%    defined to identify the candidate.                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

rng(0); % seed (reproducibility)

path = pwd;

candidate_num = 5; % identifier of the model candidate
filename = strcat('results/results_optimisation_stage1_candidate',num2str(candidate_num)); % name of the file where the results will be stored
figname = strcat('figures/figure_optimisation_stage1_candidate',num2str(candidate_num)); % name of the figure

%% loading experimental data
global ind
cd ..
load data/expDATA
A = expDATA;
ind = [2 3];
cd(path)

global mu sigma days xref

mu = A.mean;    
sigma = A.std;
days = A.days;

%% Parameters
% initial parameter values - if a parameter is to be neglected in this
% model candidate, the initial value should be set to zero
tau0 = 42.8;
tauA = 88.3; 
e = 1.1e-6; %f = 0;
m = 0.00001; %m = 0; 
d = 0.0001; %d = 0;

xref = [tau0,tauA,e,m,d]';

%% optimisation
X0 = ones(length(xref),1);
LB = zeros(1,length(xref)); 
UB = [500,500,1,1,1];

stol = 1e-8;
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','display','iter','StepTolerance',stol);

FUN = @(x) objfun_stage1(x);

tic
[x,resnorm,fval] = lsqnonlin(FUN,X0,LB,UB,options);
toc

parameters = xref.*x;

%% results
% parameter loading
param.tau0 = parameters(1); 
param.tauA = parameters(2); 
param.e = parameters(3); 
param.m = parameters(4);
param.d = parameters(5);

% simulation
opts = odeset('RelTol',1e-4,'AbsTol',1e-5);

Tf = days(end)*24;
dt = 1;
tspan = 0:dt:Tf;
Y0 = [mu(1,1); 0];

[tv,yv] = ode45(@(t,Y) odefun_stage1(t,Y,param),tspan,Y0,opts);
tsol = tv;
ysol = yv;
cells = ysol(:,1)+ ysol(:,2);

usim = interp1(tsol,cells,days*24);

N = length(mu(:,1))-1;
n = nnz(xref);

N1 = length(mu(:,1))-1;
n1 = nnz(xref);
nu1 = (N1-n1);
CHI2 = resnorm; 
AIC = CHI2 + 2*n1 + (2*n1*(n1+1))/(N1-n1-1);
BIC = CHI2 + n1*log(N1);


%% Figure
s = 15;
lw = 1.5;

fig = figure('units','normalized','outerposition',[0 0.1 0.5 0.8],'color','white');
plot(tsol,cells,'linewidth',lw,'color',[0, 0.75, 0.75]); hold on;
plot(tsol,ysol(:,1),'linewidth',lw,'color',[0 0.5 0]);
plot(tsol,ysol(:,2),'linewidth',lw,'color',[1 0 0]);
errorbar(days*24,mu(:,1),sigma(:,1),'linewidth',lw,'color',[0, 0.4470, 0.7410]);
legend('Sim Total','Sim A','Sim D','Exp','fontsize',s,'interpreter','latex','location','northwest')
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
xlabel('Time [h]','fontsize',s,'interpreter','latex')
ylabel('Number of cells','fontsize',s,'interpreter','latex')
aux = ylim;
text(350,0.92*aux(2),sprintf('%s%0.3f%s','$\chi^2 = ',CHI2,'$'),'fontsize',s,'interpreter','latex');
text(350,0.87*aux(2),sprintf('%s%0.3f%s','$AIC = ',AIC,'$'),'fontsize',s,'interpreter','latex');
text(350,0.82*aux(2),sprintf('%s%0.3f%s','$BIC = ',BIC,'$'),'fontsize',s,'interpreter','latex')

cd .. 
print(fig,figname,'-dpng');
save(filename,'parameters','xref','AIC','BIC','CHI2'); 
cd(path)


