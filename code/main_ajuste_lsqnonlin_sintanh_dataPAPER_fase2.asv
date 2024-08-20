clear all
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

filename1 = 'results_modelselection/results_f1_5.mat';

%% FASE 2
mask = ones(8,12);
mask(:,8) = 0;
mask(2,12) = 0;
mask(3,10:11) = 0;
mask(4,9) = 0;
mask(5,10:12) = 0;
mask(6,[9,12]) = 0;
mask(7,9:11) = 0;
mask(8,9:12) = 0;


%parameters
L = 25; 
n = 2;
kappa = 1e-3; 
dth = 0.95; 
DeltaD = 0.3; 
tauR = 1e9; 
lambda = 4.5e-4; 
Tmin_ini = 0;
H = 1.04;
dthA = 1400; 
DeltaDA = 10; 
Beta = 1; 

global xref2
xref2_base(1) = L; 
xref2_base(2) = n; 
xref2_base(3) = kappa; 
xref2_base(4) = dth; 
xref2_base(5) = DeltaD; 
xref2_base(6) = tauR;
xref2_base(7) = lambda;
xref2_base(8) = Tmin_ini;
xref2_base(9) = H;
xref2_base(10) = dthA;
xref2_base(11) = DeltaDA;
xref2_base(12) = Beta;

% pruebas mayo
xref5_base(1) = 160; 
xref5_base(2) = 1.7; 
xref5_base(3) = 1.3e-3; 
xref5_base(4) = 0.28; 
xref5_base(5) = 0.40; 
xref5_base(6) = 5.8e9;
xref5_base(7) = 8e-4;
xref5_base(8) = 0;
xref5_base(9) = 1.04;
xref5_base(10) = dthA;
xref5_base(11) = DeltaDA;
xref5_base(12) = Beta;

vectorAIC = zeros(1,size(mask,1));
vectorBIC = zeros(1,size(mask,1));
vectorCHI2 = zeros(1,size(mask,1));
vectorCHI2norm = zeros(1,size(mask,1));
%%
for IND = 1:8
    intento = num2str(IND);
    if IND == 8
        xref2 = xref5_base.*mask(IND,:);
    else
        xref2 = xref5_base.*mask(IND,:);
    end

    LB = zeros(1,length(xref2));
    UB = inf*ones(1,length(xref2));
    UB(12) = 1;
    UB(1) = 1.2;
    X0 = ones(1,length(xref2));
    X0(12) = 0.5;

    stol = 1e-10;
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','display','iter','StepTolerance',stol);

    FUN = @(x) objfun_f2(x,filename1);

    [x,resnorm,fval] = lsqnonlin(FUN,X0,LB,UB,options);

    parameters = xref2.*x;

    filename2 = strcat('results_modelselection/resultsBIC_f2_',intento);
    save(filename2,'parameters')

    % SIM F2
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
    Lambda = parameters(7); 
    Tmin_ini = parameters(8);
    H = parameters(9);
    dthA = parameters(10);
    DeltaDA = parameters(11);
    Beta = parameters(12);

    % model functions
    Lag = @(x) x^n/(x^n+L^n);

    Tmin = @(x) Tmin_ini + H.*x;

    Theta = @(x,y) kappa*(x>Tmin(y));
    Phi = @(x,y) (x>Tmin(y)); 

    k = log(2)/1.8; 

    A = (tau0-tauR)/tauR;
    B = (tau0+tauR)/tauR;
    if parameters(10)==0
       PIgr = @(x,y,t) 1/2*(B+A*tanh((x*Lag(t)-dth)/DeltaD));  
    else 
        PIgr = @(x,y,t) 1/2*(B+A*tanh((x*Lag(t)-dth)/DeltaD))*1/2*(1-tanh((y-dthA)/DeltaDA)); 
    end

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
            [tv,yv] = ode45(@(t,Y) odefun_f2(t,Y,param,populations(j)),tspan,Y0,opts);
            tsol = [tsol; tv(1:(end-1))];
            ysol = [ysol; yv(1:(end-1),:)];
        end 

        %add last point
        tsol = [tsol; tv(end)];
        ysol = [ysol; yv(end,:)];

        YSOL(:,:,j) = ysol;
    end

    N2 = length(populations)*(length(mu(:,ind(j)))-1);
    n2 = nnz(xref2);
    nu2 = (N2-n2);
    CHI2_normal2 = resnorm;
    AIC = CHI2_normal2 + 2*n2 + (2*n2*(n2+1))/(N2-n2-1);
    vectorAIC(IND) = AIC;
    BIC = CHI2_normal2 + n2*log(N2);
    vectorBIC(IND) = BIC;
    vectorCHI2(IND) = resnorm;

    % figure
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
        text(350,0.92*aux(2),sprintf('%s%0.3f%s','$\chi^2 = ',CHI2_normal2,'$'),'fontsize',s,'interpreter','latex')
        text(350,0.87*aux(2),sprintf('%s%0.3f%s','$\bar{\chi}^2 = ',resnorm,'$'),'fontsize',s,'interpreter','latex')
        text(350,0.82*aux(2),sprintf('%s%0.3f%s','$AIC = ',AIC,'$'),'fontsize',s,'interpreter','latex')
        text(350,0.77*aux(2),sprintf('%s%0.3f%s','$BIC = ',BIC,'$'),'fontsize',s,'interpreter','latex')
        tit = strcat('Intento ',num2str(IND));
        title(tit,'interpreter','latex','fontsize',s);
        % figname = strcat('Figures/TMZ_pruebaPEN_lsqnonlin_f1_',num2str(i),'_',intento);
        % print(fig,figname,'-dpng');
    end 
    save(filename2,'parameters','xref','resnorm','AIC','BIC');
end
save results_fase2_BIC vectorAIC vectorBIC vectorCHI2 xref2_base xref5_base