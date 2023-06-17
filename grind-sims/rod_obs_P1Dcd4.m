%% Observers for multi-model observer simulations
%

addpath("../process-observers")


%% Process and disturbance models

% This script defines the following models
%  - Gd : Process transfer function Gd(z)
%  - HDd : Input disturbance transfer function (RODD)
%  - Gpd : combined system transfer function (1 input, 1 output)
%  - Gpss : state space model of combined system

% Sample time (hrs)
Ts = 3/60;

% Continuous time process model from system identification - P1Dcd4
%              Kp                                                   
%   G(s) = ---------- * exp(-Td*s)                                  
%           1+Tp1*s                                                 
%                                                                   
%         Kp = -32.948 +/- 1.9624                                   
%        Tp1 = 0.21951 +/- 0.048178                                 
%         Td = 0.2                                                  
%
model_name = 'P1Dcd4';
s = tf('s');
Gc = -32.95 * exp(-0.2*s) / (1 + 0.2195*s);
Gc.TimeUnit = 'hours';
Gd = c2d(Gc, Ts, 'ZOH');
Gdss = ss(absorbDelay(Gd));


% RODD step disturbance process
ThetaD = 1;
PhiD = 1;
d = 1;
HDd = rodd_tf(ThetaD, PhiD, d, Ts);
HDd.TimeUnit = 'hours';

% Combined system transfer functions
Gpd = series(Gd, HDd);

% State space representation
Gpss2 = minreal(absorbDelay(ss(Gpd)));

% Discrete time state space model
A = [  0.7963        2        0        0        0        0
            0        0        1        0        0        0
            0        0        0        1        0        0
            0        0        0        0        1        0
            0        0        0        0        0        1
            0        0        0        0        0        1];
B = [0;
     0;
     0;
     0;
     0;
     1];
C = [ -3.3561       0       0       0       0       0];
D = 0;
% Alternative
% A = [       0          1;
%       -0.8899     1.8899];
% B = [0;
%      1];
% C = [-3.771 0];
% D = 0;
Gpss = ss(A, B, C, D, Ts, 'TimeUnit', 'hours');

% Dimensions
n = size(A, 1);
ny = size(C, 1);

% Check all versions are almost identical
wp = zeros(15, 1);
wp([2 6]) = [1 -1];
t_test_sim = Ts*(0:14)';
y1 = lsim(Gpd, wp, t_test_sim);
y2 = lsim(Gpss, wp, t_test_sim);
y3 = lsim(Gpss2, wp, t_test_sim);
assert(max(abs(y1 - y2)) < 0.01)
assert(max(abs(y1 - y3)) < 0.01)

% Designate which input and output variables are
% measured
u_meas = false;
y_meas = true;

% Default initial condition
x0 = zeros(n, 1);


%% Parameters for random inputs to simulations

% RODD random variable parameters
epsilon = 0.01;
sigma_wp = [0.002717 0.2717];  % [0.01 1] .* step magnitude

% Process noise standard deviation
sigma_W = [0; 0];

% Measurement noise standard deviation
sigma_M = 5;

% To check observer with no noise disturbances
%sigma_W = [0; 0];
%sigma_M = 0.001;

% Initial state of disturbance process
p0 = 0;


%% Define observers

% Observer model (without disturbance noise input)
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% Check observability of system
Qobs = obsv(A, C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Common parameters for all observers
P0 = diag([0.1*ones(1, n-1) 0.01]);  % initial error covariance
R = sigma_M^2;

% process noises for state x1
q00 = 0.1^2;

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrix
Q = diag([q00*ones(1, n-1) sigma_wp(1)^2]);
KF1 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF1');

% Kalman filter 2 - manually tuned
% Covariance matrix
Q = diag([q00*ones(1, n-1) 0.027^2]);
KF2 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF2');

% Kalman filter 3 - tuned to sigma_wp(2)
% Covariance matrix
Q = diag([q00*ones(1, n-1) sigma_wp(2)^2]);
KF3 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF3');

% Scheduled Kalman filter
% Use this to test performence of multi-model filters
Q0 = diag([q00*ones(1, n-1) 0]);
SKF = kalman_filter(A,Bu,C,Du,Ts,P0,Q0,R,'SKF');
SKF.Q0 = Q0;
SKF.Bw = Bw;
SKF.sigma_wp = sigma_wp;

% Multiple model Kalman filter
Q0 = diag([q00*ones(1, n-1) 0]);
f = 100;  % sequence history length
n_filt = 20;  % number of filters
n_min = 18;  % minimum life of cloned filters
MMKF = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,'MMKF');
