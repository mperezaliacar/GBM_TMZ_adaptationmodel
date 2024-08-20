clear 
close all
clc

tstart = tic;

N = 50; % number of sampling points per parameter
pert = linspace(-0.8,0.8,N); % vary parameters in a range of +-50% their optimal value
PERT = allcomb(pert,pert,pert,pert);

%% load experimental data
cd ..
cd ..
load expDATA/expDATA_paper_N
A = DATA_paperN;
ind = [2 3];
cd ajuste/Nivel1

mu = A.mean;    
sigma = A.std;
days = A.days;


%% load optimal (ref) value of parameters
filename1 = 'results_modelselection/results_lsqnonlin_f1_paperTea5';
filename2 = 'results_modelselection/resultsBIC_lsqnonlin_f2_UBbeta2_paperTea5';

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
num = nnz(parameters);

load(filename1,'parameters');
tau0 = parameters(1); 
tauA = parameters(2); 
e = parameters(3);
m = parameters(4); 
d = parameters(5); 

%% initialise chi2
CHI2 = zeros(length(PERT),1);

%% loop
for I=1:length(CHI2)
    disp(I)
    tic
    % perturb parameters
    kappa_pert = kappa*(1+PERT(I,1));
    dth_pert = dth*(1+PERT(I,2));
    DeltaD_pert = DeltaD*(1+PERT(I,3));
    H_pert = H*(1+PERT(I,4));
    
    % model functions
    Lag = @(x) x^n/(x^n+L^n);

    Tmin = @(x) Tmin_ini + H_pert.*x;

    Theta = @(x,y) kappa_pert*(x>Tmin(y));
    Phi = @(x,y) (x>Tmin(y)); 

    k = log(2)/1.8; %h (TMZ decay)

    A = (tau0-tauR)/tauR;
    B = (tau0+tauR)/tauR;
    PIgr = @(x,y,t) 1/2*(B+A*tanh((x*Lag(t)-dth_pert)/DeltaD_pert));  

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
            [tv,yv] = ode45(@(t,Y) odefun_f2(t,Y,param,populations(j)),tspan,Y0,opts);
            tsol = [tsol; tv(1:(end-1))];
            ysol = [ysol; yv(1:(end-1),:)];
        end 

        %add last point
        tsol = [tsol; tv(end)];
        ysol = [ysol; yv(end,:)];

        YSOL(:,:,j) = ysol;
        cells = ysol(:,1)+ysol(:,2);
        usim = interp1(tsol,cells,days*24);
        CHI2aux(j) = sum(((usim(2:end)-mu(2:end,ind(j)))./sigma(2:end,ind(j))).^2);
    end
   CHI2(I) = sum(CHI2aux);
   toc
end
CHI2matrix = reshape(CHI2,[N,N,N,N]);
CHI2matrix = permute(CHI2matrix,[4,3,2,1]);

name = strcat('loop_chi2_4param_',num2str(N),'points');
save(name,'CHI2matrix','PERT');

LOOPTIME = toc(tstart)

REFparam = [kappa,dth,DeltaD,H];
figures_contourchi2(name,REFparam,N)