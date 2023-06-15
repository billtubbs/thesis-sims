% Various observers designed for the system defined in:
%  - sys_rodin_step.m

addpath("../process-observers")

assert(exist("A", 'var'), strcat("System model not defined. ", ...
    "Run script 'sys_rodin_step.m' first."))

% Check observability of system
Qobs = obsv(A, C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Observer model without disturbance noise input
Bu = B(:, u_known);
Bw = B(:, ~u_known);
Du = D(:, u_known);

obs_model = model;
obs_model.B = Bu;
obs_model.D = Du;

% Steady-state Luenberger observer 1
% Specify poles of observer dynamics
poles = [0.8; 0.8];
LB1 = LuenbergerFilter(obs_model,poles,'LB1');

% Steady-state Luenberger observer 2
% Specify poles of observer dynamics
poles = [0.6; 0.6];
LB2 = LuenbergerFilter(obs_model,poles,'LB2');

% Specify covariance for state variable 1
% This is used by all observers
q1 = 0.01;

% Different values for covariance matrix
Q1 = diag([q1 sigma_wp{1}(1)^2]);
Q2 = diag([q1 sigma_wp{1}(2)^2]);
Q3 = diag([q1 0.1^2]);

% Covariance of output errors
R = sigma_M^2;

% Observer models for new observer functions
obs_model = model;
obs_model.B = Bu;
obs_model.D = Du;
obs_models = {obs_model, obs_model};
obs_models{1}.Q = Q1;
obs_models{1}.R = R;
obs_models{2}.Q = Q2;
obs_models{2}.R = R;

% Parameters for manually-tuned Kalman filter (KF3)
model3 = obs_models{1};  % makes copy
model3.Q = Q3;

% Steady-state Kalman filter 1 - prediction form - tuned to sigma_wp(1)
KFPSS1 = KalmanFilterPSS(obs_models{1},'KFPSS1');

% Steady-state Kalman filter 2 - prediction form - tuned to sigma_wp(2)
KFPSS2 = KalmanFilterPSS(obs_models{2},'KFPSS2');

% Steady-state Kalman filter 1 - filtering form - tuned to sigma_wp(1)
KFFSS1 = KalmanFilterFSS(obs_models{1},'KFFSS1');

% Steady-state Kalman filter 2 - filtering form  - tuned to sigma_wp(2)
KFFSS2 = KalmanFilterFSS(obs_models{2},'KFFSS2');

% Kalman filter 1 - tuned to sigma_wp(1)
P0 = 1000*eye(n);
KF1 = KalmanFilterF(obs_models{1},P0,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
P0 = 1000*eye(n);
KF2 = KalmanFilterF(obs_models{2},P0,'KF2');

% Kalman filter 3 - manually tuned
P0 = 1000*eye(n);
KF3 = KalmanFilterF(model3,P0,'KF3');

% Kalman filter 3 - prediction form - manually tuned
P0 = 1000*eye(n);
KFP3 = KalmanFilterP(model3,P0,'KFP3');

% Multiple model observer with sequence fusion based on 
% Robertson et al. (1995) paper.
label = 'MKF_SF95';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 15;  % fusion horizon
m = 1;  % maximum number of shocks
d = 3;  % spacing parameter
io = struct("u_known", u_known, "y_meas", y_meas);
MKF_SF95 = MKFObserverSF_RODD95(model,io,P0,epsilon, ...
      sigma_wp,Q0,R,f,m,d,label);

% Multiple model observer with sequence fusion #1
label = 'MKF_SF1';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
nf = 3;  % detection intervals in fusion horizon
m = 1;  % maximum number of shocks
d = 5;  % spacing parameter
io = struct("u_known", u_known, "y_meas", y_meas);
MKF_SF1 = MKFObserverSF_RODD(model,io,P0,epsilon, ...
    sigma_wp,Q0,R,nf,m,d,label);

% Multiple model observer with sequence fusion #2
label = 'MKF_SF2';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
nf = 5;  % detection intervals in fusion horizon
m = 2;  % maximum number of shocks
d = 3;  % spacing parameter
io = struct("u_known", u_known, "y_meas", y_meas);
MKF_SF2 = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nf,m,d,label);

% Multiple model observer with sequence pruning #1
label = 'MKF_SP1';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
nh = 10;  % number of hypotheses
n_min = 7;  % minimum life of cloned hypotheses
io = struct("u_known", u_known, "y_meas", y_meas);
MKF_SP1 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R, ...
    nh,n_min,label);

% Multiple model observer with sequence pruning #2
label = 'MKF_SP2';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
nh = 25;  % number of hypotheses
n_min = 21;  % minimum life of cloned hypotheses
io = struct("u_known", u_known, "y_meas", y_meas);
MKF_SP2 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R, ...
    nh,n_min,label);

observers = {LB1, LB2, KFPSS1, KFPSS2, KF1, KF2, KF3, MKF_SF95, ...
    MKF_SF1, MKF_SF2, MKF_SP1, MKF_SP2};
