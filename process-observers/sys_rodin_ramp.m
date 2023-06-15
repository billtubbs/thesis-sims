%% System model definition
%
% Discrete-time transfer function polynomial models for a
% first-order process with a RODD ramp disturbance at
% the input to the process.
% 
% Usage: run this script from a main script to define
%  - Gd : Process transfer function Gd(z)
%  - HNd : ARIMA noise disturbance transfer function HNd(z)
%  - HDd : Input disturbance RODD transfer function
%  - Gpd : combined system transfer function (3 inputs, 1 output)
%  - Gpss : state space model of combined system
%
% This script is run by run_obs_sim_spec.m during
% automated simulations.
%


%% Discrete transfer function polynomial models

% Sample time
Ts = 0.5;

% Process
omega0 = 0.3;
delta1 = 0.7;
Omega = [0 omega0];  % delay of 1
Delta = [1 -delta1];
Gd = tf(Omega, Delta, Ts);

% ARIMA noise process
% thetaN0 = 1;
% phiN1 = 0.2;
% ThetaN = [0 thetaN0];  % direct transmission
% PhiN = [1 -phiN1];
% HNd = tf(ThetaN, conv(PhiN, [1 -1]), Ts);

% RODD step disturbance process
% ThetaD = 1;
% PhiD = 1;
% d = 1;

% RODD ramp disturbance process
ThetaD = 0.001;
PhiD = 1;
d = 2;

% RODD exponential disturbance process
% ThetaD = 1;
% PhiD = [1 -0.5];
% d = 1;

HDd = rodd_tf(ThetaD, PhiD, d, Ts);

% Combined system transfer functions
Gpd = [Gd series(HDd, Gd)];

% Discrete time state space model
Gpss1 = minreal(ss(Gpd));
% A = Gpss.A;
% B = Gpss.B;
% C = Gpss.C;
% D = Gpss.D;

% Construct system manually
% Discrete time state space model
A = [0.7 1   0;
     0   1   0.001
     0   0   1];
B = [1 0;
     0 0
     0 1];
C = [0.3 0 0];
D = zeros(1, 2);
Gpss = ss(A,B,C,D,Ts);

% Check both systems are identical
nT_sim = 10;
t_sim = Ts*(0:nT_sim)';
U_sim = randn(nT_sim + 1, 2);
[Y_sim1,t_sim1,X_sim1] = lsim(Gpss1,U_sim,t_sim);
[Y_sim,t_sim,X_sim] = lsim(Gpss,U_sim,t_sim);
% TODO: Re-instate
%assert(max(abs(Y_sim1 - Y_sim)) < 1e-10)

% Designate which input and output variables are measured
u_known = [true; false];
y_meas = true;

% Dimensions
n = size(A, 1);
nu = sum(u_known);
nw = sum(~u_known);
ny = size(C, 1);

% Model parameter struct used by observers
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.D = D;
model.Ts = Ts;

% Default initial condition
x0 = zeros(n, 1);

% Parameters for random inputs
% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_wp = {[0.01 1], [0.01 1]};

% Process noise standard deviation
sigma_W = [0; 0; 0];

% Measurement noise standard deviation
sigma_M = 0.1;

% To test observer with no noise disturbances
% sigma_W = zeros(n, 1);
% sigma_M = zeros(ny, 1);

% Initial state of disturbance process
p0 = [0 0]';
