function of = objfun_stage2(x,filename1)

global mu sigma days xref2 ind 

parameters = xref2.*x;

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

cd ..
load(filename1,'parameters');
tau0 = parameters(1);
tauA = parameters(2);
e = parameters(3);
m = parameters(4);
d = parameters(5);
cd code

% model functions
Lag = @(x) x^n/(x^n+L^n);

Tmin = @(x) Tmin_ini + H.*x;

Theta = @(x,y) kappa*(x>Tmin(y));
Phi = @(x,y) (x>Tmin(y)); 

k = log(2)/1.8; %h (TMZ decay)

A = (tau0-tauR)/tauR;
B = (tau0+tauR)/tauR;
if xref2(10)==0
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

%% simulation 
populations = [2; 3];
opts = odeset('RelTol',1e-4,'AbsTol',1e-5);
T = days(end)*24;
dt = 1;

%tmz dose
hday = 24; %days-hours
daysTMZ = [1:5 28:32]-1; %days in which tmz is administered
treatment = daysTMZ*hday;

dose = 100;
CHI2aux = zeros(length(populations),length(mu)-1);

for j=1:length(populations)
    tsol = [];
    ysol = [];
    population = populations(j);
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
        [tv,yv] = ode45(@(t,Y) odefun_stage2(t,Y,param,population),tspan,Y0,opts);
        tsol = [tsol; tv(1:(end-1))];
        ysol = [ysol; yv(1:(end-1),:)];
    end 

    %add last point
    tsol = [tsol; tv(end)];
    ysol = [ysol; yv(end,:)];
    cells = ysol(:,1)+ysol(:,2);

    usim = interp1(tsol,cells,days*24);

    CHI2aux(j,:) = (usim(2:end)-mu(2:end,ind(j)))./sigma(2:end,ind(j));

end

CHI2 = [CHI2aux(1,:),CHI2aux(2,:)];
    
of = CHI2;