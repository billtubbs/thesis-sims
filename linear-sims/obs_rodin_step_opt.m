% Observers with parameters chosen by search
% See script
%  - gen_sim_specs_sim1_MKF_SP_popt.m
%  - run_obs_sims.m
% for the system defined in:
%  - sys_rodin_step.m
%

addpath("~/process-observers/")

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
obs_model.A = A;
obs_model.B = Bu;
obs_model.C = C;
obs_model.D = Du;
obs_model.Ts = Ts;

% Specify covariance for state variable 1
% This is used by all observers
q1 = 0.01;

% Different values for covariance matrix
Q1 = diag([q1 sigma_wp{1}(1)^2]);
Q2 = diag([q1 sigma_wp{1}(2)^2]);
Q3 = diag([q1 0.1^2]);

% Covariance of output errors
R = sigma_M^2;

% Observer models for multiple-model observers
obs_models = {obs_model, obs_model};
obs_models{1}.Q = Q1;
obs_models{1}.R = R;
obs_models{2}.Q = Q2;
obs_models{2}.R = R;

% Parameters for manually-tuned Kalman filter (KF3)
model3 = obs_models{1};  % makes copy
model3.Q = Q3;

% Kalman filter 1 - tuned to sigma_wp(1)
P0 = 1000*eye(n);
KF1 = KalmanFilterF(obs_models{1},P0,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
P0 = 1000*eye(n);
KF2 = KalmanFilterF(obs_models{2},P0,'KF2');

% Kalman filter 3 - manually tuned
P0 = 1000*eye(n);
KF3 = KalmanFilterF(model3,P0,'KF3');

% Multiple model observer with sequence fusion based on 
% Robertson et al. (1995) paper.
label = 'MKF_SF95';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 5;  % fusion horizon
m = 1;  % maximum number of shocks
d = 1;  % spacing parameter
io.u_known = u_known;
io.y_meas = y_meas;
MKF_SF95 = MKFObserverSF_RODD95(model,io,P0,epsilon, ...
      sigma_wp,Q0,R,f,m,d,label);

% Multiple model observer with sequence fusion based on 
% Robertson et al. (1995) paper.
label = 'MKF_SF95_2';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 10;  % fusion horizon
m = 2;  % maximum number of shocks
d = 1;  % spacing parameter
io.u_known = u_known;
io.y_meas = y_meas;
MKF_SF95_2 = MKFObserverSF_RODD95(model,io,P0,epsilon, ...
      sigma_wp,Q0,R,f,m,d,label);

% Multiple model observer with sequence fusion based on 
% Robertson et al. (1998) paper.
label = 'MKF_SF1';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
nf = 3;  % number of detection intervals in fusion horizon
m = 1;  % maximum number of shocks
d = 2;  % spacing parameter
io.u_known = u_known;
io.y_meas = y_meas;
MKF_SF1 = MKFObserverSF_RODD(model,io,P0,epsilon, ...
      sigma_wp,Q0,R,nf,m,d,label);

% Multiple model observer with sequence pruning #1
label = 'MKF_SP1';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
nh = 10;  % number of hypotheses
n_min = 7;  % minimum life of cloned hypotheses
io.u_known = u_known;
io.y_meas = true(ny, 1);
MKF_SP1 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R, ...
    nh,n_min,label);

% Multiple model observer with sequence pruning #2
label = 'MKF_SP2';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
nh = 25;  % number of hypotheses
n_min = 21;  % minimum life of cloned hypotheses
MKF_SP2 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R, ...
    nh,n_min,label);

observers = {KF1, KF2, KF3, MKF_SF95, MKF_SF95_2, MKF_SF1, MKF_SP1, MKF_SP2};
