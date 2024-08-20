%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       sensitivity_analysis.m                            %
%                                                                         %
%    This script carries out a local sensitivity analysis around the      %
%    values specified in parameters.m                                     %
%                                                                         %
%    It returns the barplot of the sensitivity index of the parameters    %
%    as well as a results file with the value of the indexes              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc

%% load experimental data
cd ..
load data/expDATA
A = DATA_paperN;
ind = [2 3];
cd code

mu = A.mean;    
sigma = A.std;
days = A.days;

%%
cd ..
load results/results_model_simulation
YSOL_ref = YSOL;
cells_ref = YSOL_ref(:,1,2:3)+ YSOL_ref(:,2,2:3);
uref = cells_ref;
cd code

%% load value of parameters
parameters;

% ref_values
L_ref = L;
n_ref = n;
kappa_ref = kappa;
dth_ref = dth;
DeltaD_ref = DeltaD;
tauR_ref = tauR;
Lambda_ref = Lambda;
H_ref = H;
dthA_ref = dthA;
DeltaDA_ref = DeltaDA;
Beta_ref = Beta;

num_PoI = 8; % parameters of interest (L,n,kappa,dth,DeltaD,tauR,Lambda,H)
ind_PoI = [1,2,3,4,5,6,7,9];
delta = 0.1; % perturb 10% of value
PoI_mask = delta*eye(8);

%%
EE = zeros(num_PoI,length(tsol),3);
EE_norm = zeros(num_PoI,3);

for I=1:num_PoI
    % perturb parameters if necessary
    L = L_ref*(1+PoI_mask(1,I));
    n = n_ref*(1+PoI_mask(2,I));  
    kappa = kappa_ref*(1+PoI_mask(3,I));
    dth = dth_ref*(1+PoI_mask(4,I));
    DeltaD = DeltaD_ref*(1+PoI_mask(5,I)); 
    tauR = tauR_ref*(1+PoI_mask(6,I)); 
    Lambda = Lambda_ref*(1+PoI_mask(7,I)); 
    H = H_ref*(1+PoI_mask(8,I));

    % model functions
    Lag = @(x) x^n/(x^n+L^n);

    Tmin = @(x) Tmin_ini + H.*x;

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

    %%%%% simulation
    opts = odeset('RelTol',1e-4,'AbsTol',1e-5);

    T = days(end)*24;
    dt = 1;

    %tmz dose
    hday = 24; %days-hours
    daysTMZ = [1:5 28:32]-1; %days in which tmz is administered
    treatment = daysTMZ*hday;
    dose = 100;
    populations = [2 3];

    CHI2aux = zeros(length(populations),1);

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
        cells = ysol(:,1)+ysol(:,2);
        usim = cells;
        
        % calculate elementary effects (EE)
        EE(I,:,j) = (usim-uref(:,:,j))/(delta);
        EE_norm(I,j) = EE(I,end,j)/usim(end);
    end
    EE(I,:,3) = 1/2*(EE(I,:,1)+EE(I,:,2)); 
    EE_norm(I,3) = 1/2*(EE_norm(I,1)+EE_norm(I,2));
end
filename = strcat('results/local_sens_analysis_delta',num2str(delta*100));
cd ..
save(filename,'EE','delta');
cd code

%% figure
names = {'$L$','$n$','$\kappa$','$S_{th}$','$\Delta S$','$\tau_\mathrm{s}$','$\lambda$','$\gamma$'};
titles = {'Sensibles','Resistentes','Promedio'};
colors = {[0.4940, 0.1840, 0.5560],[0.9290, 0.6940, 0.1250],[0.4660, 0.6740, 0.1880]};
paramcolors = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],[0.75, 0.75, 0],[0, 0.5, 0],[0.75, 0, 0.75]};	

posX = 3;
posY = 3;
width = 12;
heigth = 10;
lw = 1.5;
s = 10;


for i=1:3
    fig = figure();
    set(gcf,'units','centimeters','position',[posX,posY,width,heigth]);
    bar(EE_norm(:,i),'FaceColor',colors{i})
    xticks(1:8);
    title(titles{i},'interpreter','latex','fontsize',s)
    xticklabels(names)
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15);
    ylim([-1.2,1.2]);
    figname = strcat('figures/barplot',titles{i});
    cd ..
    print(fig,figname,'-dpng','-r600');
    cd code
end
