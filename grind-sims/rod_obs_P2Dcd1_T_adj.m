%% Observers for multi-model observer simulations
%
% This script contains an adjustable model which can be 
% specified with different parameters.
%
% Before running this script define the following variables:
% Kp : default -35.94
% Tp1 : default 0.2350
% Tp2 : default 0.1608
% thetap : default 0.05
% epsilon : default 0.01
% step_mag : default 0.2717
% b : default 100
% where sigma_wp = [1/b 1] .* step_mag
%
% See original model script for details:
%  - rod_obs_P2Dcd1_T.m

addpath("../process-observers")


%% Process and disturbance models

% This script defines the following models
%  - Gd : Process transfer function Gd(z)
%  - HDd : Input disturbance transfer function (RODD)
%  - Gpd : combined system transfer function (1 input, 1 output)
%  - Gpss : state space model of combined system

% Sample time (hrs)
Ts = 3/60;

% Continuous time process model from system identification 
% with true outputs - P2Dcd1_T
%
%                 Kp
%   G(s) = ----------------- * exp(-Td*s)
%          (1+Tp1*s)(1+Tp2*s)
% 
%          Kp = -35.935 +/- 0.11198
%         Tp1 = 0.23503 +/- 0.017791
%         Tp2 = 0.16076 +/- 0.01594
%          Td = 0.05
%
model_name = 'P2Dcd1_T_adj';
s = tf('s');
Gc = Kp * exp(-thetap * s) / ((1 + Tp1*s) * (1 + Tp2*s));
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
Gpss = minreal(absorbDelay(ss(Gpd)));
model = struct();
model.A = Gpss.A;
model.B = Gpss.B;
model.C = Gpss.C;
model.D = Gpss.D;
model.Ts = Gpss.Ts;

% Discrete time state space model (P1D_c4)
% A = [ 2.5411  -1.0667   0.5923        1
%            2        0        0        0
%            0      0.5        0        0
%            0        0        0        0];
% B = [   0;
%         0;
%         0;
%         1];
% C = [      0 -0.50028 -0.84026        0];
% D = 0;
% Gpss = ss(A, B, C, D, Ts, 'TimeUnit', 'hours');

% Check structure matches requirements for simulations
assert(isequal(model.B, [0 0 0 1]'))

% Dimensions
[n, ~, ny] = check_dimensions(model.A, model.B, model.C, model.D);

% Check all versions are almost identical
wp = zeros(15, 1);
wp([2 6]) = [1 -1];
t_test_sim = Ts*(0:14)';
y1 = lsim(Gpd, wp, t_test_sim);
y2 = lsim(Gpss, wp, t_test_sim);
%y3 = lsim(Gpss2, wp, t_test_sim);
assert(max(abs(y1 - y2)) < 0.01)
%assert(max(abs(y1 - y3)) < 0.01)

% Designate which input and output variables are
% measured
u_known = false;
y_meas = true;

% Default initial condition
x0 = zeros(n, 1);


%% Parameters for random inputs to simulations

% RODD random variable parameters
%epsilon = 0.01;
sigma_wp = {[sigma_wp_1 sigma_wp_1*b]};

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

% Define which inputs and outputs are measured/known
io = struct();
io.u_known = u_known;
io.y_meas = y_meas;

% Observer model (without disturbance noise input)
obs_model = model;  % makes a copy
obs_model.B = model.B(:, u_known);
obs_model.D = model.D(:, u_known);
obs_model.nu = sum(u_known);

% Disturbance input (used by SKF observer)
Bw = model.B(:, ~u_known);
nu = sum(u_known);
nw = sum(~u_known);

% Check observability of system
Qobs = obsv(model.A, model.C);
unobs = length(model.A) - rank(Qobs);
assert(unobs == 0);

% Common parameters for all observers
P0 = diag([0.1*ones(1, n-1) 0.01]);  % initial error covariance
R = sigma_M^2;

% process noises for state x1
q00 = 0.01^2;

% Kalman filter 1 - tuned to sigma_wp(1)
obs_model1 = obs_model;
obs_model1.Q = diag([q00*ones(1, n-1) sigma_wp{1}(1)^2]);
obs_model1.R = R;
KF1 = KalmanFilterF(obs_model1,P0,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
obs_model2 = obs_model;
obs_model2.Q = diag([q00*ones(1, n-1) sigma_wp{1}(2)^2]);
obs_model2.R = R;
KF2 = KalmanFilterF(obs_model2,P0,'KF2');

% Kalman filter 3 - manually tuned
obs_model3 = obs_model;
obs_model3.Q = diag([q00*ones(1, n-1) 0.027^2]);
obs_model3.R = R;
KF3 = KalmanFilterF(obs_model3,P0,'KF3');

% Switching models
obs_models = {obs_model1, obs_model2};

% Scheduled Kalman filter
% Use this to test performance of multi-model filters
SKF = SKFObserver(obs_models,P0,'SKF');

% Multiple model observer - sequence fusion
Q0 = diag([q00*ones(1, n-1) 0]);
f = 20;  % fusion horizon
m = 2;  % maximum number of shocks
d = 5;  % spacing parameter
MKF_SF95 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,'MKF_SF95');

% Multiple model observer - sequence fusion
Q0 = diag([q00*ones(1, n-1) 0]);
nf = 4;  % detection intervals in fusion horizon
m = 2;  % maximum number of shocks
d = 12;  % spacing parameter
MKF_SF1 = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nf,m,d,'MKF_SF1');

% Multiple model observer - sequence pruning
Q0 = diag([q00*ones(1, n-1) 0]);
nh = 20;  % number of filters
n_min = 18;  % minimum life of cloned filters
MKF_SP1 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,'MKF_SP');

observers = {KF1, KF2, KF3, SKF, MKF_SF95, MKF_SF1, MKF_SP1};
