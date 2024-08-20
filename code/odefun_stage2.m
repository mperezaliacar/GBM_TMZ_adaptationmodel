function f = odefun_stage2(t,U,p,pop)

% load parameters
tau0 = p.tau0;
tauA = p.tauA;
e = p.e;
m = p.m;
d = p.d;

Beta = p.Beta;

if pop == 2
    Lambda = 0;
elseif pop == 3
    Lambda = p.Lambda;
end

PIgr = p.PIgr;
Theta = p.Theta;
Phi = p.Phi;

k = p.k;

N0 = 1000;

Fgr = N0^(1/3)*(U(1)+U(2))^(2/3);

% normoxic cells
f(1,1) = 1/tau0*PIgr(U(4),U(5),t)*Fgr-1/tauA*U(1)-e*U(1)*U(2)-m*U(1);

% dead cells
f(2,1) = 1/tauA*U(1)-d*U(2);

% TMZ
f(3,1) = -k*U(3);

% damage
f(4,1) = Theta(U(3),U(5))*U(3)*Phi(U(3),U(5))-Beta*1/(U(1)*tau0)*PIgr(U(4),U(5),t)*Fgr*U(4)-Lambda*U(4);

% damage AC
f(5,1) = U(4);

end 