% Simulate observers in Simulink with S-functions
%
% Simulates the same system and set of observers as in
% in test_MKF_example.m and checks the results are
% identical
%
% NOTE: This test does not work if run as a test script.
% To run tests involving Simulink models requires the 
% sltest package which I do not have. Therefore, run 
% this like a regular Matlab file (not using runtests)
%

clear all

%show_plots = true;
show_plots = false;

% See this Simulink model file:
sim_model = 'MKF_example_sim';

% Generate randomly-occurring shocks
% Reset random number generator
seed = 22;
rng(seed)

% Load system model
% First order SISO system with one input disturbance
sys_rodin_step

% Sequence length
nT = 100;

% Generate random shock signal
[Wp, alpha] = sample_random_shocks(nT+1, epsilon, sigma_wp{1}(2), ...
    sigma_wp{1}(1));

% Other inputs to system
X0 = zeros(n,1);
t = Ts*(0:nT)';
U = zeros(nT+1,1);
U(t>=5) = 1;
V = sigma_M*randn(nT+1, 1);

% Simulate system in MATLAB
[Y,T,X] = lsim(Gpss,[U Wp],t,X0);
Ym = Y + V;  % measurement

if show_plots
    figure(2)
    subplot(2,1,1)
    plot(t,Y,t,Ym); grid on
    ylabel('y(k) and y_m(k)')
    legend('y(k)', 'ym(k)')
    title('Outputs')
    subplot(2,1,2)
    stairs(t, [U Wp]); grid on
    xlabel('k')
    ylabel('u(k) and wp(k)')
    legend('u(k)', 'wp(k)')
    title('Inputs')
end

% Inputs to simulink model
inputs.U = [t U];
inputs.V = [t V];
inputs.Wp = [t Wp];

% Observer model without disturbance noise input
Bu = B(:, u_known);
Bw = B(:, ~u_known);
Du = D(:, u_known);

% Kalman filter parameters
Q = diag([0.01^2 0.1^2]);
R = 0.1^2;
obs_model = model;
obs_model.B = Bu;
obs_model.D = Du;
obs_model.Q = Q;
obs_model.R = R;

% Prediction form
KFPSS = KalmanFilterPSS(obs_model,'KFPSS');

% Filtering form
KFFSS = KalmanFilterFSS(obs_model,'KFFSS');
% Copy - used to test MKFObserver_sfunc
KFFSS2 = KalmanFilterFSS(obs_model,'KFFSS2');

% Dynamic Kalman filters
P0 = eye(n);

% Prediction form
KFP = KalmanFilterP(obs_model,P0,'KFP');

% Filtering form
KFF = KalmanFilterF(obs_model,P0,'KFF');

% Copy - used to test MKFObserver_sfunc
KFF2 = KalmanFilterF(obs_model,P0,'KFF2');

% Switching Kalman filter - scheduled
obs_models = {obs_model, obs_model};
P0 = eye(n);
obs_models{1}.Q = diag([0.01^2 sigma_wp{1}(1)^2]);
obs_models{2}.Q = diag([0.01^2 sigma_wp{1}(2)^2]);
seq = alpha' + 1;
SKF = SKFObserverS(obs_models,P0,seq,"SKF");

% Multi-model observer 1 - sequence fusion
P0 = eye(n);
Q0 = diag([0.01^2 0]);
R = 0.1^2;
nf = 5;  % fusion horizon
m = 1;  % maximum number of shocks
d = 3;  % spacing parameter
io.u_known = u_known;
io.y_meas = y_meas;
MKF1 = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp,Q0,R, ...
    nf,m,d,'MKF1');

% Multi-model observer 2 - sequence pruning
P0 = eye(n);
Q0 = diag([0.01^2 0]);
R = 0.1^2;
nh = 5;  % number of filters
n_min = 2;  % minimum life of cloned filters
io.u_known = u_known;
io.y_meas = true(ny, 1);
MKF2 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,'MKF2');

if show_plots
    open(sim_model)
end
fprintf("Running Simulink simulation...\n")
%sim_out = sim(sim_model, 'StopTime', string(nT*Ts), ...
%    'ReturnWorkspaceOutputs', 'on');
sim_out = sim(sim_model, 'StopTime', string(nT*Ts));

% Check built-in Simulink Kalman filter block estimates are same as
% steady-state observer object - prediction form
assert(max(abs(sim_out.X_hat_KF.Data - sim_out.X_hat_KFPSS.Data), [], [1 2]) < 5e-15)

% Check all filter estimates are close to true system states
assert(sqrt(mean((sim_out.X_hat_KFPSS.Data - sim_out.X.Data).^2, [1 2])) < 0.30)
assert(sqrt(mean((sim_out.X_hat_KFP.Data - sim_out.X.Data).^2, [1 2])) < 0.31)
assert(sqrt(mean((sim_out.X_hat_KFFSS.Data - sim_out.X.Data).^2, [1 2])) < 0.23)
assert(sqrt(mean((sim_out.X_hat_KFF.Data - sim_out.X.Data).^2, [1 2])) < 0.23)
assert(sqrt(mean((sim_out.X_hat_SKF.Data - sim_out.X.Data).^2, [1 2])) < 0.15)
%assert(sqrt(mean((sim_out.X_hat_MKF1.Data - sim_out.X.Data).^2, [1 2])) < 0.99)
assert(sqrt(mean((sim_out.X_hat_MKF2.Data - sim_out.X.Data).^2, [1 2])) < 0.18)

% Check KalmanFilterF_sfunc.m and MKFObserver_sfunc.m
% produce identical results with KFF, KFFSS
assert(isequal(sim_out.X_hat_KFFSS.Data, sim_out.X_hat_KFFSS2.Data))
assert(isequal(sim_out.X_hat_KFF.Data, sim_out.X_hat_KFF2.Data))

% Load simulation results produced by test_MKF_example.m
filename = 'MKF_example.csv';
results_dir = 'results';
sim_results = readtable(fullfile(results_dir, filename));

% This loads simulation results for the following observers:
% - Xk_est_KFPSS, Yk_est_KFPSS
% - Xk_est_KFFSS, Yk_est_KFFSS
% - Xk_est_KFP, Yk_est_KFP
% - Xk_est_KFF, Yk_est_KFF
% - Xk_est_SKF, Yk_est_SKF
% - Xk_est_MKF1, Yk_est_MKF1
% - Xk_est_MKF2, Yk_est_MKF2

% Check all Simulink observer estimates are same as simulation results
% from file

% Steady-state Kalman filter - prediction form
assert(max(abs(sim_out.X_hat_KFPSS.Data - ...
    sim_results{:, {'Xk_est_KFPSS_1', 'Xk_est_KFPSS_2'}}), [], [1 2]) < 1e-12)

% Steady-state Kalman filter - filtering form
assert(max(abs(sim_out.X_hat_KFFSS.Data - ...
    sim_results{:, {'Xk_est_KFFSS_1', 'Xk_est_KFFSS_2'}}), [], [1 2]) < 1e-12)

% Dynamic Kalman filter - prediction form
assert(max(abs(sim_out.X_hat_KFP.Data - ...
    sim_results{:, {'Xk_est_KFP_1', 'Xk_est_KFP_2'}}), [], [1 2]) < 1e-12)

% Dynamic Kalman filter - filtering form
assert(max(abs(sim_out.X_hat_KFF.Data - ...
    sim_results{:, {'Xk_est_KFF_1', 'Xk_est_KFF_2'}}), [], [1 2]) < 1e-12)

% Switching Kalman filter
assert(max(abs(sim_out.X_hat_SKF.Data - ...
    sim_results{:, {'Xk_est_SKF_1', 'Xk_est_SKF_2'}}), [], [1 2]) < 1e-12)

% Multiple-model RODD observer - SF
assert(max(abs(sim_out.X_hat_MKF1.Data - sim_results{:, {'Xk_est_MKF1_1', 'Xk_est_MKF1_2'}}), [], [1 2]) < 1e-12)

% Multiple-model RODD observer - SP
assert(max(abs(sim_out.X_hat_MKF2.Data - sim_results{:, {'Xk_est_MKF2_1', 'Xk_est_MKF2_2'}}), [], [1 2]) < 1e-12)

disp("Simulations complete")
