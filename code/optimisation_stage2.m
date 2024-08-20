%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         optimisation_stage2.m                           %
%                                                                         %
%    This script performs the STAGE 2 of the optimisation of the model    %
%    parameters, that is, those related to temozolomide response, for     %
%    a given model candidate (defined by setting to zero the value of     %
%    some parameters) and computes the resulting values of the chi^2,     %
%    BIC and AIC that can be used to perform model selection.             %
%                                                                         %
%    The user must specify the model candidate as well as give a number   %
%    for each candidate (candidate_num) to identify the results.          %       
%                                                                         %
%    The results are stored in the results folder as:                     %
%    "results_optimisation_stage2_candidateX.mat" with X the number       %
%    defined to identify the candidate.                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

rng(0); % seed (reproducibility)

path = pwd;

candidate_num = 5; % identifier of the model candidate
filename2 = strcat('results/results_optimisation_stage2_candidate',num2str(candidate_num)); % name of the file where the results will be stored
figname = strcat('figures/figure_optimisation_stage2_candidate',num2str(candidate_num)); % name of the figure
population_names = {'Sensitive','Resistant'};

%%%%%%%%%%%%%%%%%%% Exp data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ind
cd ..
load data/expDATA
A = DATA_paperN;
ind = [2 3];
cd(path)

global mu sigma days xref

mu = A.mean;    
sigma = A.std;
days = A.days;

optimal_candidate_stage1 = 5;
filename1 = strcat('results/results_optimisation_stage1_candidate',num2str(candidate_num)); % name of the file where the optimal results of stage 1 were stored

%% Parameters
% initial parameter values - if a parameter is to be neglected in this
% model candidate, the initial value should be set to zero (except for dth
% A, that should be set to infinity
L = 160; 
n = 1.7;
kappa = 1.3e-3; 
dth = 0.28; 
DeltaD = 0.4; 
tauR = 5.8e9; 
lambda = 8e-4; 
Tmin_ini = 0;
H = 1.04;
dthA = 1400; dthA = inf;
DeltaDA = 10; DeltaDA = 0;
Beta = 1;  Beta = 0;

global xref2
xref2(1) = L; 
xref2(2) = n; 
xref2(3) = kappa; 
xref2(4) = dth; 
xref2(5) = DeltaD; 
xref2(6) = tauR;
xref2(7) = lambda;
xref2(8) = Tmin_ini;
xref2(9) = H;
xref2(10) = dthA;
xref2(11) = DeltaDA;
xref2(12) = Beta; 


%% optimisation
LB = zeros(1,length(xref2));
UB = inf*ones(1,length(xref2));
UB(1) = 1.2;
UB(12) = 1;
X0 = ones(1,length(xref2));
X0(12) = 0.5;

stol = 1e-10;
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','display','iter','StepTolerance',stol);

FUN = @(x) objfun_stage2(x,filename1);

[x,resnorm,fval] = lsqnonlin(FUN,X0,LB,UB,options);

parameters2 = xref2.*x;

%% results
% parameter loading
L = parameters2(1);
n = parameters2(2); 
kappa = parameters2(3); 
dth = parameters2(4);
DeltaD = parameters2(5); 
tauR = parameters2(6);
Lambda = parameters2(7); 
Tmin_ini = parameters2(8);
H = parameters2(9);
dthA = parameters2(10);
DeltaDA = parameters2(11);
Beta = parameters2(12);

% recover parameters of stage 1
cd ..
load(filename1,'parameters')
tau0 = parameters(1); 
tauA = parameters(2); 
e = parameters(3); 
m = parameters(4);
d = parameters(5);
cd(path)

% model functions
Lag = @(x) x^n/(x^n+L^n);

Tmin = @(x) Tmin_ini + H.*x;

Theta = @(x,y) kappa*(x>Tmin(y));
Phi = @(x,y) (x>Tmin(y)); 

k = log(2)/1.8; 

A = (tau0-tauR)/tauR;
B = (tau0+tauR)/tauR;
PIgr = @(x,y,t) 1/2*(B+A*tanh((x*Lag(t)-dth)/DeltaD))*1/2*(1-tanh((y-dthA)/DeltaDA)); 

%loading
param.Lag = Lag;
param.Theta = Theta;
param.Phi = Phi;
param.Tmin = Tmin;
param.k = k;
param.PIgr = PIgr;
param.Beta = Beta;
param.Lambda = Lambda;

param.tau0 = tau0;
param.tauA = tauA;
param.e = e;
param.m = m;
param.d = d;

% simulation
opts = odeset('RelTol',1e-4,'AbsTol',1e-5);

T = days(end)*24;
dt = 1;

%tmz dose
hday = 24; 
daysTMZ = [1:5 28:32]-1; %days in which tmz is administered
treatment = daysTMZ*hday;

dose = 100;
populations = [2 3];

for j=1:length(populations)
    tsol = [];
    ysol = [];
    for i=1:length(treatment)
        % temporal bounds
        Ti = treatment(i);
        if i == length(treatment)
            Tf = T;
        else 
            Tf = treatment(i+1);
        end
        tspan = Ti:dt:Tf;
        % initial conditions
        Y0 = zeros(5,1);
        if i == 1
            Y0(1) = mu(1,3); %cells
            Y0(2) = 0;
            Y0(3) = dose;
            Y0(4) = 0; %damage
            Y0(5) = 0; %damage ac
        else 
            Y0(1) = yv(end,1);
            Y0(2) = yv(end,2);
            Y0(3) = dose;
            Y0(4) = yv(end,4);
            Y0(5) = yv(end,5);
        end 
        [tv,yv] = ode45(@(t,Y) odefun_stage2(t,Y,param,populations(j)),tspan,Y0,opts);
        tsol = [tsol; tv(1:(end-1))];
        ysol = [ysol; yv(1:(end-1),:)];
    end 

    %add last point
    tsol = [tsol; tv(end)];
    ysol = [ysol; yv(end,:)];

    YSOL(:,:,j) = ysol;
end

N2 = length(populations)*(length(mu(:,ind(j)))-1);
n2 = nnz(xref2)-sum(isinf(xref2));
nu2 = (N2-n2);
CHI2 = resnorm;
AIC = CHI2 + 2*n2 + (2*n2*(n2+1))/(N2-n2-1);
BIC = CHI2 + n2*log(N2);


%% Figure
s = 15;
lw = 1.5;
colors = {[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250]};
for i=1:length(populations)
    fig = figure('units','normalized','outerposition',[0.5 0.1 0.5 0.8],'color','white');
    cells = YSOL(:,1,i)+YSOL(:,2,i);
    plot(tsol,cells,'linewidth',lw,'color',[0, 0.75, 0.75]); hold on;
    plot(tsol,YSOL(:,1,i),'linewidth',lw,'color',[0 0.5 0]);
    plot(tsol,YSOL(:,2,i),'linewidth',lw,'color',[1 0 0]);
    errorbar(days*24,mu(:,ind(i)),sigma(:,ind(i)),'linewidth',lw,'color',colors{i});
    legend('Sim Total','Sim A','Sim D','Exp','fontsize',s,'interpreter','latex','location','northwest')
    ax = gca;
    ax.FontSize = floor(0.9*s);
    ax.TickLabelInterpreter = 'latex';
    xlabel('Time [h]','fontsize',s,'interpreter','latex')
    ylabel('Number of cells','fontsize',s,'interpreter','latex')
    aux = ylim;
    text(350,0.92*aux(2),sprintf('%s%0.3f%s','$\chi^2 = ',CHI2,'$'),'fontsize',s,'interpreter','latex')
    text(350,0.87*aux(2),sprintf('%s%0.3f%s','$AIC = ',AIC,'$'),'fontsize',s,'interpreter','latex')
    text(350,0.82*aux(2),sprintf('%s%0.3f%s','$BIC = ',BIC,'$'),'fontsize',s,'interpreter','latex')
    title(population_names{i},'interpreter','latex','fontsize',s);
    cd ..
    print(fig,strcat(figname,'_',population_names{i}),'-dpng');
    cd(path)
end 
cd ..
save(filename2,'parameters2','xref','CHI2','AIC','BIC');
cd(path)
