% System model definition
%
% Discrete-time switching system for testing GPB1 and GPB2
% state estimators.
% 
% Usage: run this script from a main script to define
%  - A1, B1, C1, D1
%  - A2, B2, C2, D2
%  - Ts
%  - m1, m2, models
%
% This system is used for testing in the following scripts:
%  - test_MKFObservers_JS.m 
%  - MKFObserverGPB.m
%

% Define system

% Sample period
Ts = 0.5;

% Discrete time state space models
% Model #1
A1 = 0.7;
B1 = 1;
C1 = 0.3;
D1 = 0;
%Gpss1 = ss(A1,B1,C1,D1,Ts);

% Model #2
A2 = 0.9;
B2 = 1;
C2 = -0.3;  % -ve gain!
D2 = 0;
%Gpss2 = ss(A2,B2,C2,D2,Ts);

% Two system models as structs
m1 = struct();
m1.A = A1;
m1.B = B1;
m1.C = C1;
m1.D = D1;
m1.Ts = Ts;
m2 = struct();
m2.A = A2;
m2.B = B2;
m2.C = C2;
m2.D = D2;
m2.Ts = Ts;

% Cell array containing both models
models = {m1, m2};
[nj, n, nu, ny, Ts] = check_models(models);

% Input disturbance variance
%sigma_w = 0.1;
sigma_w = 0;

% Process noise std. dev.
sigma_W = [0; 0];

% Measurement noise std. dev.
sigma_M = 0.1;  % Note: above 0.2 numerical error occurs in GPB2 algorithm

% Initial condition
x0 = 0;

