clear
close all
filename1 = 'results_modelselection/results_lsqnonlin_f1_paperTea5';
filename2 = 'results_modelselection/resultsBIC_lsqnonlin_f2_UBbeta2_paperTea5';
cd ..
cd ..
%% datos fase 1
load expDATA/expDATA_paper_N
A = DATA_paperN;

%%
cd ajuste/Nivel1
mu = A.mean;    
sigma = A.std;
days = A.days;
ind = [2 3];
lw = 1.5;
s = 15;
%% results
% SIM F2
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
dthA = parameters(10); %dthA = 2000;
DeltaDA = parameters(11); %DeltaDA = 90;
Beta = parameters(12); %Beta = 0.1;
num = nnz(parameters);

load(filename1,'parameters');
tau0 = parameters(1); 
tauA = parameters(2); 
e = parameters(3);
m = parameters(4); 
d = parameters(5); 

% model functions
Fgr = @(x,y) (x+y)^(2/3);

Lag = @(x) x^n/(x^n+L^n);

Tmin = @(x) Tmin_ini + H.*x;

Theta = @(x,y) kappa*(x>Tmin(y));
Phi = @(x,y) (x>Tmin(y)); 

k = log(2)/1.8; %h (TMZ decay)

A = (tau0-tauR)/tauR;
B = (tau0+tauR)/tauR;
if dthA==0
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
param.Fgr = Fgr;
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

% poblacion control
opts = odeset('RelTol',1e-4,'AbsTol',1e-5);

dt = 1;
tspan = 0:dt:T;
Y0 = [mu(1,1); 0];

[tv,yv] = ode45(@(t,Y) odefun_f1(t,Y,param),tspan,Y0,opts);
tsolcontrol = tv;
ysolcontrol = yv;

% poblaciones tratadas
dose = 100;
populations = [2 3];
tic;
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
save results_optimal YSOL tsol
toc;
%% figure
CHI2aux = zeros(length(populations),1);
colors = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250]};


% poblaciones tratadas
for i=1:length(populations)
    fig = figure('units','normalized','outerposition',[0.5 0.1 0.25 0.5],'color','white');
    cells = YSOL(:,1,i)+YSOL(:,2,i);
    plot(tsol,cells,'linewidth',lw,'color',[0, 0.75, 0.75]); hold on;
    plot(tsol,YSOL(:,1,i),'linewidth',lw,'color',[0 0.5 0]);
    plot(tsol,YSOL(:,2,i),'linewidth',lw,'color',[1 0 0]);
    errorbar(days*24,mu(:,ind(i)),sigma(:,ind(i)),'linewidth',lw,'color',colors{ind(i)});
%     legend('Sim Total','Sim A','Sim D','Exp','fontsize',s,'interpreter','latex','location','northwest')
    ax = gca;
    ax.FontSize = floor(0.9*s);
    ax.TickLabelInterpreter = 'latex';
    xlabel('Time [h]','fontsize',s,'interpreter','latex')                                                     
    ylabel('Number of cells','fontsize',s,'interpreter','latex')
    aux = ylim;
    usim = interp1(tsol,cells,days*24);
    N = 2*(length(mu(:,1))-1);
    CHI2aux(i) = sum(((usim(2:end)-mu(2:end,ind(i)))./sigma(2:end,ind(i))).^2);
end 
chi2 = sum(CHI2aux);
text(350,0.92*aux(2),sprintf('%s%0.3f%s','$\chi^2 = ',chi2,'$'),'fontsize',s,'interpreter','latex')
N2 = N;
n2 = num;
nu2 = (N2-n2);
CHI2_normal2 = chi2;
AIC = CHI2_normal2 + 2*n2 + (2*n2*(n2+1))/(N2-n2-1)
BIC = CHI2_normal2 + n2*log(N2)

% Tmin = @(x) Tmin_ini + H*x;
% 
% fig = figure('units','normalized','outerposition',[0.1 0.1 0.5 0.5],'color','white');
% subplot(1,3,1)
% plot(tsol,YSOL(:,3,2)); hold on;
% plot(tsol,Tmin(YSOL(:,5,2)));
% subplot(1,3,2)
% plot(tsol,YSOL(:,4,2)); hold on;
% yline(dth,'k--','linewidth',lw,'DisplayName','Threshold')
% subplot(1,3,3)
% plot(tsol,YSOL(:,5,2)); hold on;
% yline(dthA,'k--','linewidth',lw,'DisplayName','Threshold')