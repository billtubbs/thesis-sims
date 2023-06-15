% Test function run_simulation_obs.m

clear all

% Specify path to observer functions
addpath("../process-observers")

% Data file used in test
results_dir = "../process-observers/results";

% SISO system example from GEL-7029 course, homework 12.
% See file /gel-7029/homework/hw12/hw12_p3_kalman.m

sys_test_siso2


%% Define and simulate Kalman filters

% Observer parameters
W = 0.5; % estimate of Wp used in the filter design
V = 0.8; % estimate of Vp used in the filter design
P0 = eye(n)*1000; % Initialize covariance matrix
model.Q = diag(repmat(W,n,1));
model.R = diag(repmat(V,ny,1));

% Define dynamic Kalman filter using kalman_filter function
KF = KalmanFilterP(model,P0,"KF");

assert(all(isnan(KF.K)))
assert(isequal(KF.Pkp1, P0))
assert(isequal(KF.xkp1_est, zeros(2, 1)))
assert(KF.ykp1_est == 0)

observers = {KF};

% number of points to simulate
nT = 100;

% seed random number generator
rng(0)

% Measurement noise for the whole simulation
V = (sqrt(Rp)*randn(1, nT))';

% Process noise for the whole simulation
W = (sqrt(Qp)*randn(2, nT))';

% Intialize system (at k = 0)
x0 = zeros(n, 1);
x = x0;

% Input signal - pseudo-random binary sequence
warnId = 'Controllib:estimation:initialPRBSSequence';
warnStruct = warning('off',warnId);
U = idinput(nT, 'PRBS', [0 0.5]);
warning(warnStruct);

% Designate measured/unmeasured inputs and outputs variables
u_meas = true;
y_meas = true;
alpha = randn(nT, 1);  % not used in this test
Pd = nan(nT, 1);  % not used in this test

% Define state transition function and measurement
% function of system (used by run_simulation_obs)
[f, h, params] = make_trans_funcs_ssd(A, B, C, D);

input_data = table(alpha, U, V, W, Pd);

% Simulate system and observers
sim_out = run_simulation_obs(Ts, input_data, f, h, params, x0, ...
    u_meas, observers);

% Display results
%head(sim_out.data)

% Verify results by comparing with Kalman_Filter.mlx

filename = 'hw12_p3_kalman_sim_benchmark.csv';

warnId = 'MATLAB:table:ModifiedAndSavedVarnames';
warnStruct = warning('off',warnId);
bench_sim_results = readtable(fullfile(results_dir, filename));
warning(warnStruct);

% Display benchmark simulation results
%head(bench_sim_results)

% Check simulation inputs are identical
assert(isequal( ...
    round(sim_out.data{1:10, {'U', 'V', 'W'}}, 6), ...
    round(bench_sim_results{1:10, {'u_t_', 'v_t_', 'w_t__1', 'w_t__2'}}, 6) ...
))

% Check states
assert(isequal( ...
    round(sim_out.data{1:10, 'X'}, 6), ...
    round(bench_sim_results{1:10, {'X_t__1', 'X_t__2'}}, 6) ...
))

% Check state estimates
assert(isequal( ...
    round(sim_out.data{1:10, 'X_est'}, 6), ...
    round(bench_sim_results{1:10, {'X_e_t__1', 'X_e_t__2'}}, 6) ...
))

% plot results

% t = sim_out.data{:, 't'};
% U = sim_out.data{:, 'U'};
% X = sim_out.data{:, 'X'};
% Y_m = sim_out.data{:, 'Y_m'};
% X_est = sim_out.data{:, 'X_est'};
% Y_est = sim_out.data{:, 'Y_est'};
% 
% figure(1); clf
% 
% subplot(411);
% plot(t, Y_m, 'k', t, Y_est, 'r', 'Linewidth', 2)
% legend('Process output', 'KF estimates')
% ylabel('$y$')
% grid on
% 
% subplot(412);
% plot(t, X(:,1), 'k', t, X_est(:,1), 'r', 'Linewidth', 2)
% legend('Process state', 'KF estimates')
% ylabel('$x_1$')
% grid on
% 
% subplot(413);
% plot(t, X(:,2), 'k', t, X_est(:,2), 'r', 'Linewidth', 2);
% legend('Process state', 'KF estimates')
% ylabel('$x_2$')
% grid on
% 
% subplot(414);
% stairs(t', U', 'Linewidth', 2);
% xlabel('Time [s]');
% ylabel('$u_1$')
% grid on

