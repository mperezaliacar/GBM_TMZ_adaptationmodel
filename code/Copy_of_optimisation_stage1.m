%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         optimisation_stage1.m                           %
%                                                                         %
%    This script provides the result of the model simulation with the     %
%    parameter values specified in the file parameters.m, returning       % 
%    figures with both experimental (contained in the file expDATA.m)     %
%    and simulated data for the three populations.                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

rng(0); % seed (reproducibility)

%%%%%%%%%%%%%%%%%%% Exp data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ind
cd ..
cd ..
load expDATA/expDATA_paper_N
A = DATA_paperN;
ind = [2 3];
cd ajuste/Nivel1

global mu sigma days xref

mu = A.mean;    
sigma = A.std;
days = A.days;

%% define model candidates
mask = ones(8,5);
mask(2,5) = 0;
mask(3,4) = 0;
mask(4,3) = 0;
mask(5,4:5) = 0;
mask(6,[3,5]) = 0;
mask(7,3:4) = 0;
mask(8,3:5) = 0;


%% Fase 1: control
% initial parameter values
tau0 = 42.8;
tauA = 88.3; 
e = 1.1e-6; 
m = 0.00001; 
d = 0.0001; 

xref_base = [tau0,tauA,e,m,d];

vectorAIC = zeros(1,size(mask,1));
vectorBIC = zeros(1,size(mask,1));
vectorCHI2 = zeros(1,size(mask,1));

for IND = 1:size(mask,1)
    intento = num2str(IND);
    xref = (xref_base.*mask(IND,:))';
    
    X0 = ones(length(xref_base),1);
    LB = zeros(1,length(xref_base)); 
    UB = [500,500,1,1,1];
    tol = 1e-2;
    stol = 1e-8;
    
    FUN = @(x) objfun_f1(x);

    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','display','iter','StepTolerance',stol);

    tic
    [x,resnorm,fval] = lsqnonlin(FUN,X0,LB,UB,options);
    toc
    
    parameters = xref.*x;

    filename1 = strcat('results_modelselection/results_f1_',intento);

    %% results
    % SIM F1
    %loading
    param.tau0 = parameters(1); 
    param.tauA = parameters(2); 
    param.e = parameters(3); 
    param.m = parameters(4);
    param.d = parameters(5);

    opts = odeset('RelTol',1e-4,'AbsTol',1e-5);

    Tf = days(end)*24;
    dt = 1;
    tspan = 0:dt:Tf;
    Y0 = [mu(1,1); 0];

    [tv,yv] = ode45(@(t,Y) odefun_f1(t,Y,param),tspan,Y0,opts);
    tsol = tv;
    ysol = yv;
    cells = ysol(:,1)+ ysol(:,2);

    usim = interp1(tsol,cells,days*24);

    N = length(mu(:,1))-1;
    n = nnz(xref);
    
    N1 = length(mu(:,1))-1;
    n1 = nnz(xref);
    nu1 = (N1-n1);
    CHI2_normal = resnorm; 
    AIC = CHI2_normal + 2*n1 + (2*n1*(n1+1))/(N1-n1-1);
    BIC = CHI2_normal + n1*log(N1);
    
    vectorAIC(IND) = AIC;
    vectorBIC(IND) = BIC;
    vectorCHI2(IND) = CHI2_normal;

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
    text(350,0.92*aux(2),sprintf('%s%0.3f%s','$\chi^2 = ',CHI2_normal,'$'),'fontsize',s,'interpreter','latex');
    text(350,0.87*aux(2),sprintf('%s%0.3f%s','$\bar{\chi}^2 = ',resnorm,'$'),'fontsize',s,'interpreter','latex'); % cambio resnorm por fval
    text(350,0.82*aux(2),sprintf('%s%0.3f%s','$AIC = ',AIC,'$'),'fontsize',s,'interpreter','latex');
    text(350,0.77*aux(2),sprintf('%s%0.3f%s','$BIC = ',BIC,'$'),'fontsize',s,'interpreter','latex')

    % figname = strcat('Figures/Control_lsqnonlin_f1_',intento);
    % print(fig,figname,'-dpng');
    save(filename1,'parameters','xref','resnorm','AIC','BIC'); % cambio resnorm por fval
end 
save results_fase1_BIC_pruebalsqnonlin vectorAIC vectorBIC vectorCHI2 xref_base
