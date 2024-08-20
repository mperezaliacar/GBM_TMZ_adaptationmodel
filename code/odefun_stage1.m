function f = odefun_stage1(t,U,p)

% load parameters
tau0 = p.tau0;
tauA = p.tauA;
e = p.e;
m = p.m;
d = p.d;

N0 = 1000;

Fgr = N0^(1/3)*(U(1)+U(2))^(2/3);

% normoxic cells
f(1,1) = 1/tau0*Fgr-1/tauA*U(1)-e*U(1)*U(2)-m*U(1);

% dead cells
f(2,1) = 1/tauA*U(1)-d*U(2);

end 