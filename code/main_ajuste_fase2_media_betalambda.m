clear 
close all
clc

rng(0)

global xref_repair

%%%%%%%%%%%%%%%%%%% Exp data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..
cd ..
load expDATA/expDATA_ind
name = 'paperTea';
cd ajuste/Nivel1

filename1 = 'results_modelselection/results_lsqnonlin_f1_paperTea5.mat';
filename2 = 'results_modelselection/results_lsqnonlin_f2_mediatratados';

N = size(datos_tratados,1);
days = time;

mu = mean(datos_tratados,1);
sigma = std(datos_tratados,1);

load(filename1)
tau0 = parameters(1); 
tauA = parameters(2); 
e = parameters(3); 
m = parameters(4);
d = parameters(5);

load(filename2)
L = parameters(1);
n = parameters(2);
kappa = parameters(3);
dth = parameters(4);
DeltaD = parameters(5);
tauR = parameters(6);
H = parameters(7);

xref_repair(1) = 6.9e-4; %lambda
xref_repair(2) = 1; %beta

% model
Lag = @(x) x^n/(x^n+L^n);

Tmin = @(x) H.*x;

Theta = @(x,y) kappa*(x>Tmin(y));
Phi = @(x,y) (x>Tmin(y)); 

k = log(2)/1.8; %h (TMZ decay)

A = (tau0-tauR)/tauR;
B = (tau0+tauR)/tauR;

PIgr = @(x,y,t) 1/2*(B+A*tanh((x*Lag(t)-dth)/DeltaD));  

% loading
param.Lag = Lag;
param.Theta = Theta;
param.Phi = Phi;
param.Tmin = Tmin;
param.k = k;
param.PIgr = PIgr;

param.tau0 = tau0;
param.tauA = tauA;
param.e = e;
param.m = m;
param.d = d;

% optimisation
LB = [0 0];
UB = [inf 1];
X0 = [1 0.8];  

stol = 1e-10;
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','display','iter','StepTolerance',stol);

CHI2aux = zeros(N,1);
vectorBeta = zeros(N,1);
vectorLambda = zeros(N,1);

% figure
s = 8;
lw = 1;

fig = figure('units','normalized','outerposition',[0 0 1 1],'color','white');

for j=1:N
    uexp = datos_tratados(j,:);
    FUN_repair = @(x) objfun_betalambda(x,uexp,param);
    [x,~,fval] = lsqnonlin(FUN_repair,X0,LB,UB,options);
    CHI2aux(j) = fval;
    vectorLambda(j) = x(1)*xref_repair(1);
    vectorBeta(j) = x(2)*xref_repair(2);   
    
    %% simulation 
    load(filename1,'parameters')
    tau0 = parameters(1); 
    tauA = parameters(2); 
    e = parameters(3); 
    m = parameters(4);
    d = parameters(5);

    load(filename2,'parameters');
    L = parameters(1);
    n = parameters(2);
    kappa = parameters(3);
    dth = parameters(4);
    DeltaD = parameters(5);
    tauR = parameters(6);
    H = parameters(7);
    
    Lambda = x(1)*xref_repair(1); %Lambda = 0;
    Beta = x(2)*xref_repair(2); %Beta = parameters(9); 

    % model functions
    Lag = @(x) x^n/(x^n+L^n);

    Tmin = @(x) H.*x;

    Theta = @(x,y) kappa*(x>Tmin(y));
    Phi = @(x,y) (x>Tmin(y)); 

    k = log(2)/1.8; %h (TMZ decay)

    A = (tau0-tauR)/tauR;
    B = (tau0+tauR)/tauR;
    PIgr = @(x,y,t) 1/2*(B+A*tanh((x*Lag(t)-dth)/DeltaD));  

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

    opts = odeset('RelTol',1e-4,'AbsTol',1e-5);

    T = days(end)*24;
    dt = 1;

    %tmz dose
    hday = 24; %days-hours
    daysTMZ = [1:5 28:32]-1; %days in which tmz is administered
    treatment = daysTMZ*hday;

    dose = 100;
    populations = [2 3];

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
            Y0(1) = 1000; %cells
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
        [tv,yv] = ode45(@(t,Y) odefun_f2_media(t,Y,param),tspan,Y0,opts);
        tsol = [tsol; tv(1:(end-1))];
        ysol = [ysol; yv(1:(end-1),:)];
    end 

    %add last point
    tsol = [tsol; tv(end)];
    ysol = [ysol; yv(end,:)];

    N2 = length(uexp)-1;
    n2 = 2;
    nu2 = (N2-n2);
    CHI2_normal2 = fval*nu2;
    BIC = CHI2_normal2 + n2*log(N2);

    subplot(8,12,j)
    cells = ysol(:,1) + ysol(:,2);
    plot(tsol,cells,'linewidth',lw,'color',[0, 0.75, 0.75]); hold on;
    plot(tsol,ysol(:,1),'linewidth',lw,'color',[0 0.5 0]);
    plot(tsol,ysol(:,2),'linewidth',lw,'color',[1 0 0]);
    plot(days*24,uexp,'linewidth',lw,'color',[0.8500, 0.3250, 0.0980]);
%     legend('Sim Total','Sim A','Sim D','Exp','fontsize',s,'interpreter','latex','location','northwest')
    ax = gca;
    ax.FontSize = floor(0.9*s);
    ax.TickLabelInterpreter = 'latex';
%     xlabel('Time [h]','fontsize',s,'interpreter','latex')
%     ylabel('Number of cells','fontsize',s,'interpreter','latex')
    aux = ylim;
%     text(350,0.92*aux(2),sprintf('%s%0.3f%s','$\chi^2 = ',CHI2_normal2,'$'),'fontsize',s,'interpreter','latex')
%     text(350,0.87*aux(2),sprintf('%s%0.3f%s','$\bar{\chi}^2 = ',fval,'$'),'fontsize',s,'interpreter','latex')
%     text(350,0.82*aux(2),sprintf('%s%0.3f%s','$BIC = ',BIC,'$'),'fontsize',s,'interpreter','latex')

end

plot(vectorLambda(1:69),vectorBeta(1:69),'*r'); hold on; plot(vectorLambda(70:end),vectorBeta(70:end),'*b'); xlabel('Lambda'); ylabel('Beta');

save data_media_betalambda_1 vectorLambda vectorBeta
