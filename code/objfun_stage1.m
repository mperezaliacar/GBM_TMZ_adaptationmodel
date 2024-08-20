function of = objfun_stage1(x)

global mu sigma days xref 

parameters = xref.*x;

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

[tv,yv] = ode45(@(t,Y) odefun_stage1(t,Y,param),tspan,Y0,opts);
tsol = tv;
ysol = yv;

cells = ysol(:,1)+ysol(:,2);

usim = interp1(tsol,cells,days*24);

N = length(mu(:,1))-1;
n = nnz(xref);

CHI2 = ((usim(2:end)-mu(2:end,1))./sigma(2:end,1));

of = CHI2;