%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            parameters.m                                 %
%                                                                         %
%    This script contains the value of the model parameters, to be        %
%    defined by the user for running model_simulation.m                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parameters spheroid growth
tau0 = 42.5; % growth characteristic time [h]
tauA = 83.49; % death characteristic time [h] 
e = 1.1e-6; % growth inhibition coefficient [cell^-1·h^-1]
m = 0; % shedding rate [h^-1]
d = 0; % dead cells disappearance rate [h^-1]

%% parameters TMZ response
L = 192; % lag time [h]
n = 1.4; % Lag sharpness coefficient [-] 
kappa = 1.4e-3; % stress acquisition coefficient [h^-1]
dth = 0.29; % stress threshold [-]
DeltaD = 0.37; % stress spread parameter [-]
tauR = 5.8e9; % cytostatic growth characteristic time [h]
Lambda = 8.4e-4; % stress decay rate [h^-1]
Tmin_ini = 0; % initial TMZ threshold [muM]
H = 1; % adaptation coefficient [muM·cell^-1]
dthA = inf; % accumulated stress threshold [h]
DeltaDA = 1; % accumulates stress spread parameter [h]
Beta = 0; % inheritance coefficient [-]
