% SISO system example from GEL-7029 homework 12, Problem 3
% See file /gel-7029/homework/hw12/hw12_p3_kalman.m

% Sampling period
Ts = 2;

% Plant
A = [1.5 -0.6; 1 0];
B = [1; 0];
C = [-0.05 0.5];
D = 0;

% Model parameter struct used by observers
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.D = D;
model.Ts = Ts;

% Dimensions
n = size(A, 1); % number of states
nu = size(B, 2);  % number of inputs
ny = size(C, 1); % number of outputs

% Covariance of process noise
Qp = [0.5 0; 0 0.5];

% Variance of measurement noise
Rp = 0.8;

% To check observer with no noise disturbances
%Qp = zeros(2, 2); % process noise
%Rp = 0; % measurement noise