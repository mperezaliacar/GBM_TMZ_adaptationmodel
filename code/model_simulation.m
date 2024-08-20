%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         model_simulation.m                              %
%                                                                         %
%    This script provides the result of the model simulation with the     %
%    parameter values specified in the file parameters.m, returning       % 
%    figures (in the figures folder) with both experimental (contained    % 
%    in the file expDATA.m) and simulated data for the three              %
%    populations:                                                         %
%       - Results_control.png: control population                         %
%       - Results_TMZsens.png: TMZ-sensitive population                   %
%       - Results_TMZres.png: TMZ-resistant population                    %
%                                                                         %
%    It also returns a .mat file (results_model_simulation.mat) in the    %
%    results forlder with the cell, species and internal variables        %
%    values along the simulation time                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

%% load experimental data
cd ..
load data/expDATA
A = DATA_paperN;
x = A.days;
y = A.mean;
err = A.std;
days = A.days;
IND = [1 2 3]; 
cd code

%% load parameters
parameters

%% model 
% model functions
Fgr = @(x,y) (x+y)^(2/3);

Lag = @(x) x^n/(x^n+L^n);

Tmin = @(x) Tmin_ini + H.*x;

Theta = @(x,y) kappa*(x>Tmin(y));
Phi = @(x,y) (x>Tmin(y)); 

k = log(2)/1.8; %h (TMZ decay)

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
param.Fgr = Fgr;
param.Beta = Beta;
param.Lambda = Lambda;

param.tau0 = tau0;
param.tauA = tauA;
param.e = e;
param.m = m;
param.d = d;

opts = odeset('RelTol',1e-4,'AbsTol',1e-5);


%tmz dose
hday = 24; %days-hours
daysTMZ = [1:5 28:32]-1; %days in which tmz is administered
treatment = daysTMZ*hday;

T = days(end)*24;
dt = 1;
tspan = 0:dt:T;
Y0 = [y(1,1); 0];
YSOL = zeros(length(tspan),5,3);

[tv,yv] = ode45(@(t,Y) odefun_stage1(t,Y,param),tspan,Y0,opts);
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
        [tv,yv] = ode45(@(t,Y) odefun_stage2(t,Y,param,populations(j)),tspan,Y0,opts);
        tsol = [tsol; tv(1:(end-1))];
        ysol = [ysol; yv(1:(end-1),:)];
    end 

    %add last point
    tsol = [tsol; tv(end)];
    ysol = [ysol; yv(end,:)];
    
    YSOL(:,:,j+1) = ysol;
end
cd ..
cd results
save results_model_simulation YSOL tsol

%% FIGURES
cd ..
cd figures
% colors
TMZ = [0.4940, 0.1840, 0.5560];	
TMZ_BIG = [0.9290, 0.6940, 0.1250];
DMSO_BIG = [0, 0.4470, 0.7410];

colors = {DMSO_BIG,TMZ,TMZ_BIG};

colorsim = {[0,0,1],[247, 37, 202]/255,[250, 156, 32]/255};
color_marker = [0.8500, 0.3250, 0.0980];

LIMS = [14e4, 2e4, 4e4];

posX = 10;
posY = 10;
width = 5.75;
heigth = 6;
s = 9;
lw = 1;


%% Figures of comparison between simulated and experimental data
for i=1:3
    fig = figure();
    set(gcf,'units','centimeters','position',[posX,posY,width,heigth],'color','white');
    cells = YSOL(:,1,i)+YSOL(:,2,i);
    errorbar(days*24,y(:,i),err(:,i),'linewidth',lw,'color',colors{i},'DisplayName','Experimental results'); hold on;
    plot(tsol,cells,'linewidth',1.4*lw,'color',colorsim{i},'linestyle','-','DisplayName','Total number of cells'); 
    plot(tsol,YSOL(:,1,i),'linewidth',1,'color',[0.4660, 0.6740, 0.1880],'linestyle','-','DisplayName','Alive cells');
    plot(tsol,YSOL(:,2,i),'linewidth',1,'color',[1 0 0],'linestyle','-','DisplayName','Dead cells');
    xx = [0 1 2 3 4 28 29 30 31 32]*24;
    yy = -y(1,1)/2*ones(size(xx));

    leg = legend('fontsize',s,'interpreter','latex','location','northwest');
    leg.ItemTokenSize = [11,14];
    leg.Box = 'off';
    ax = gca;
    ax.FontSize = floor(0.9*s);
    ax.TickLabelInterpreter = 'latex';
    xlabel('Time [h]','interpreter','latex','fontsize',s)
    ylabel('Cell number','interpreter','latex','fontsize',s)
    xlim([0 1.01*x(end)*24]);
    ylim([0 LIMS(i)]);    

    figname = {'Results_control','Results_TMZsens','Results_TMZres'};
    print(fig,figname{i},'-dpng','-r600');
end

%% Figures of TMZ and TMZ threshold
posX = 10;
posY = 10;
width = 8;
heigth = 7; height = 5; 

TT = Tmin(YSOL(:,5,3));
ind = find(floor(TT)==400);

fig = figure(5);
set(gcf,'units','centimeters','position',[posX,posY,width,height]);

Xfill = [tsol(1:ind)', tsol(ind), 0, 0];
Yfill = [TT(1:ind)', TT(ind), TT(ind), 0]; 

plot(tsol,TT,'linewidth',0.8*lw,'linestyle','--','color',[0.6350, 0.0780, 0.1840],'DisplayName','$T_\mathrm{min}$'); hold on;
plot(tsol,YSOL(:,3,3),'linewidth',lw,'color',[0, 0.75, 0.75],'DisplayName','$T$'); 
fill(Xfill,Yfill,[0.6350, 0.0780, 0.1840],'FaceAlpha',0.2,'EdgeColor','none','DisplayName','Drug effect area'); 

xlim([0 1.01*x(end)*24]);  
ylim([0 350])
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
xlabel('Time [h]','interpreter','latex','fontsize',s)
ylabel('TMZ concentration [$\mu$M]','interpreter','latex','fontsize',s)
leg = legend('fontsize',s,'interpreter','latex','location','northeast');
leg.ItemTokenSize = [11,14];
leg.Box = 'on';
print(fig,'TMZ_threshold','-dpng','-r600');

cd ..
cd code