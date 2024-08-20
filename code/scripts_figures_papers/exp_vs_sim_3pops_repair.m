clear 
close all
clc

IND = 1;

%% Results DEF
filename1 = 'results_modelselection/results_lsqnonlin_f1_paperTea5';
filename2 = strcat('results_modelselection/resultsBIC_lsqnonlin_f2_repair_paperTea',num2str(IND));
cd ..
cd ..
cd ..
load expDATA/expDATA_paper_N
A = DATA_paperN;
x = A.days;
y = A.mean;
err = A.std;
days = A.days;
cd ajuste/Nivel1

%%
%% results
% SIM F2
load(filename2,'parameters');
L = parameters(1); %L = 100;
n = parameters(2); 
kappa = parameters(3); %kappa = 0.001;
dth = parameters(4); %dth = 0.8;
DeltaD = parameters(5); %DeltaD = 0.35;
tauR = parameters(6); 
Lambda = parameters(7); %Lambda = 1e-5;
Tmin_ini = parameters(8);
H = parameters(9); %H = 0.08;
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
Y0 = [y(1,1); 0];
YSOL = zeros(length(tspan),5,3);

[tv,yv] = ode45(@(t,Y) odefun_f1(t,Y,param),tspan,Y0,opts);
tsolcontrol = tv;
ysolcontrol = yv;

YSOL(:,1:2,1) = ysolcontrol;

% poblaciones tratadas
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
            Y0(1) = y(1,3); %cells
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
        [tv,yv] = ode45(@(t,Y) odefun_f2_repair(t,Y,param,populations(j),IND),tspan,Y0,opts);
        tsol = [tsol; tv(1:(end-1))];
        ysol = [ysol; yv(1:(end-1),:)];
    end 

    %add last point
    tsol = [tsol; tv(end)];
    ysol = [ysol; yv(end,:)];
    
    YSOL(:,:,j+1) = ysol;
end


%% FIGURES
% colors
TMZ = [0.4940, 0.1840, 0.5560];	
TMZ_BIG = [0.9290, 0.6940, 0.1250];
DMSO_BIG = [0, 0.4470, 0.7410];

colors = {DMSO_BIG,TMZ,TMZ_BIG};
% colors{3} = [0, 0.75, 0.75];

colorsim = {[0, 0, 1],[235,39,242]/255,[219,124,0]/255};
colorsim = {[0,0,1],[247, 37, 202]/255,[250, 156, 32]/255};
% colorsim{3} = [0, 0.55, 0.55];

color_marker = [0.8500, 0.3250, 0.0980];

LIMS = [14e4, 2e4, 4e4];

posX = 10;
posY = 10;
width = 5.75;
heigth = 6;
s = 9;
lw = 1;

cd scripts_figures_papers

%% Figura - ajuste de las 3 poblaciones juntas en inglés
% width = 8;
% heigth = 7;
% fig = figure(4);
% set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
% figname = strcat('figure_repair_intento_',num2str(IND));
% namexp = {'Exp: Control','Exp: TMZ-S','Exp: TMZ-R'};
% namsim = {'Sim: Control','Sim: TMZ-S','Sim: TMZ-R'};
% for i=1:3
%     errorbar(days*24,y(:,i),err(:,i),'linewidth',lw,'color',colors{i},'DisplayName',namexp{i}); hold on;
% end
% for i=1:3
%     cells = YSOL(:,1,i)+YSOL(:,2,i);
%     plot(tsol,cells,'linewidth',1.4*lw,'color',colorsim{i},'linestyle','-','DisplayName',namsim{i}); 
% end
% 
% leg = legend('fontsize',s,'interpreter','latex','location','northwest','numcolumns',2);
% leg.ItemTokenSize = [11,14];
% leg.Box = 'off';
% ax = gca;
% ax.FontSize = floor(0.9*s);
% ax.TickLabelInterpreter = 'latex';
% xlabel('Time [h]','interpreter','latex','fontsize',s)
% ylabel('Cell number','interpreter','latex','fontsize',s)
% xlim([0 1.01*x(end)*24]);    
% 
% print(fig,figname,'-dpng','-r600');

%% Figura - graphical abstract
% width = 7;
% heigth = 7;
% fig = figure(4);
% set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
% figname = strcat('figure_repair_intento_',num2str(IND));
% namexp = {'Exp: Control','Exp: TMZ-S','Exp: TMZ-R'};
% namsim = {'Sim: Control','Sim: TMZ-S','Sim: TMZ-R'};
% for i=1:3
%     errorbar(days*24,y(:,i),err(:,i),'linewidth',lw,'color',colors{i},'DisplayName',namexp{i}); hold on;
% end
% for i=1:3
%     cells = YSOL(:,1,i)+YSOL(:,2,i);
%     plot(tsol,cells,'linewidth',1.4*lw,'color',colorsim{i},'linestyle','-','DisplayName',namsim{i}); 
% end
% 
% % leg = legend('fontsize',s,'interpreter','latex','location','northwest','numcolumns',2);
% % leg.ItemTokenSize = [11,14];
% % leg.Box = 'off';
% ax = gca;
% ax.FontSize = floor(0.9*s);
% ax.TickLabelInterpreter = 'latex';
% xlabel('Time [h]','interpreter','latex','fontsize',s)
% ylabel('Cell number','interpreter','latex','fontsize',s)
% xlim([0 1.01*x(end)*24]);    
% ylim([0 13e4])
% 
% print(fig,'graphical_abstract_CIBM','-dpng','-r600');

%% Figura - abstract IUTAM JAJ
width = 7;
heigth = 4;
fig = figure(4);
set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
figname = strcat('figure_abstrac_IUTAM_JAJ',num2str(IND));
namexp = {'Exp: Control','Exp: TMZ-S','Exp: TMZ-R'};
namsim = {'Sim: Control','Sim: TMZ-S','Sim: TMZ-R'};
for i=2:3
    errorbar(days*24,y(:,i),err(:,i),'linewidth',lw,'color',colors{i},'DisplayName',namexp{i}); hold on;
end
for i=2:3
    cells = YSOL(:,1,i)+YSOL(:,2,i);
    plot(tsol,cells,'linewidth',1.4*lw,'color',colorsim{i},'linestyle','-','DisplayName',namsim{i}); 
end

leg = legend('fontsize',s,'interpreter','latex','location','northwest','numcolumns',1);
leg.ItemTokenSize = [11,14];
leg.Box = 'off';
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
xlabel('Time [h]','interpreter','latex','fontsize',s)
ylabel('Cell number','interpreter','latex','fontsize',s)
xlim([0 1.01*x(end)*24]);    

print(fig,figname,'-dpng','-r600');
