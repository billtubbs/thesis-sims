% Test the following multi-model observer classes:
%  - MKFObserverAMM, MKFObserverGPB1, and MKFObserverGPB2

clear

seed = 0;
rng(seed)

addpath("../plot-utils")

data_dir = 'data';


%% Simulation test - SISO system

% Load switching system definition
sys_js2_siso

% Default measurement noise std. dev. (see sys_switch_JS.m)
% sigma_M = 0.1;
% TODO: above 0.2 numerical error occurs in GPB2 algorithm

% Simulation settings
nT = 60;
t = Ts*(0:nT)';

% Inputs
%U = (idinput(size(t)) + 1)/2;
U = zeros(nT+1,1);
U(t>2) = 1;
V = sigma_M * randn(size(t));

% Actual system switching sequence
%Gamma = int8(rand(nT+1, 1) > T(1, 1));
Gamma = int8(zeros(nT+1, 1));
Gamma(t>=10, 1) = 1;

% Simulate switching system
[X, Y, Ym] = run_simulation_sys(models,U,V,Gamma,nT);

% % Plot of system inputs and outputs
% figure(1); clf
% make_iodplot(Y,Ym,t,[U Gamma+1],{'$u(k)$','$r(k)$'},{'$y(k)$','$y_m(k)$'})

% Observer parameters
P0 = 10000;
m1.Q = 0.01;
m1.R = 0.1^2;
m2.Q = m1.Q;  % same for both modes
m2.R = m1.R;  % same for both modes
models = {m1, m2};
nj = numel(models);

% Initial condition
x0 = 0.5;

% Standard Kalman filters
KF1 = KalmanFilterF(m1,P0,'KF1',x0);
KF2 = KalmanFilterF(m2,P0,'KF2',x0);

% Define scheduled MKF filter
seq = Gamma' + 1;
SKF = SKFObserverS(models,P0,seq,"SKF_S");

assert(strcmp(SKF.type, "SKF_S"))
assert(isequal(SKF.models, models))
assert(isequal(SKF.P0, P0))
assert(isequaln(SKF.seq, seq))
assert(strcmp(SKF.label, "SKF_S"))
assert(isequal(SKF.x0, zeros(n, 1)))
assert(SKF.n == n)
assert(SKF.nu == nu)
assert(SKF.ny == ny)
assert(SKF.nj == nj)
assert(isequal(SKF.xkp1_est, zeros(n, 1)))
assert(isequal(SKF.Pkp1, P0))
assert(isequal(SKF.rk, seq(:, 1)))
assert(isequaln(SKF.xk_est, nan(n, 1)))
assert(isequaln(SKF.Pk, nan(n)))
assert(isequaln(SKF.yk_est, nan(ny, 1)))
assert(isequaln(SKF.Kf, nan(n, ny)))
assert(isequaln(SKF.Sk, nan(ny)))
assert(isequaln(SKF.nf, size(seq, 2)))
assert(isequaln(SKF.i, 0))
assert(isequaln(SKF.i_next, 1))

% Transition probabilities
T = [0.95 0.05; 0.01 0.99];
assert(all(sum(T, 2) == 1))

% System indicator sequences
seq1 = {
    zeros(1, nT+1);
    [zeros(1, 20) ones(1, nT+1-20)];  % equal to Gamma'
    [zeros(1, 40) ones(1, nT+1-40)];
    ones(1, nT+1);
 };
assert(isequal(seq1{2}, Gamma'))

% Define MKF observers

% 1. AMM
MKF_AMM1 = MKFObserverAMM(models,P0,"MKF_AMM1",x0);

assert(strcmp(MKF_AMM1.type, "MKF_AMM"))
assert(isequal(MKF_AMM1.models, models))
assert(isequal(MKF_AMM1.Ts, Ts))
assert(isequal(MKF_AMM1.P0, P0))
assert(strcmp(MKF_AMM1.label, "MKF_AMM1"))
assert(MKF_AMM1.nh == 2)
assert(MKF_AMM1.n == n)
assert(MKF_AMM1.nu == nu)
assert(MKF_AMM1.ny == ny)
assert(MKF_AMM1.nj == 2)
assert(isequal(MKF_AMM1.xkp1_est, x0))
assert(isequal(MKF_AMM1.Pkp1, P0))
assert(isequaln(MKF_AMM1.xk_est, nan(n, 1)))
assert(isequaln(MKF_AMM1.Pk, nan(n)))
assert(isequaln(MKF_AMM1.yk_est, nan(ny, 1)))
assert(isequaln(MKF_AMM1.x0, x0))
assert(isequal(MKF_AMM1.p_seq_g_Yk_init, [0.5 0.5]'))

% Also with initial prior probability values
p_seq_g_Yk_init = [0.6; 0.4];
MKF_AMM1x0p0 = MKFObserverAMM(models,P0,"MKF_GPB1x0",x0, ...
    p_seq_g_Yk_init);
assert(isequaln(MKF_AMM1x0p0.x0, x0))
assert(isequaln(MKF_AMM1x0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF_AMM1x0p0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))

% 2. GPB1
MKF_GPB1 = MKFObserverGPB1(models,P0,T,"MKF_GPB1",x0);

assert(strcmp(MKF_GPB1.type, "MKF_GPB1"))
assert(isequal(MKF_GPB1.models, models))
assert(isequal(MKF_GPB1.Ts, Ts))
assert(isequal(MKF_GPB1.P0, P0))
assert(isequal(MKF_GPB1.T, T))
assert(strcmp(MKF_GPB1.label, "MKF_GPB1"))
assert(MKF_GPB1.nh == 2)
assert(MKF_GPB1.n == n)
assert(MKF_GPB1.nu == nu)
assert(MKF_GPB1.ny == ny)
assert(MKF_GPB1.nj == 2)
assert(isequal(MKF_GPB1.T, T))
assert(isequal(MKF_GPB1.xkp1_est, x0))
assert(isequal(MKF_GPB1.Pkp1, P0))
assert(isequaln(MKF_GPB1.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB1.Pk, nan(n)))
assert(isequaln(MKF_GPB1.yk_est, nan(ny, 1)))
assert(isequaln(MKF_GPB1.rk, [1 2]'))
assert(isequaln(MKF_GPB1.p_yk_g_seq_Ykm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.p_seq_g_Ykm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.p_rk_g_rkm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.p_seq_g_Ykm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.x0, x0))
assert(isequal(MKF_GPB1.p_seq_g_Yk_init, [0.5 0.5]'))
assert(isequaln(MKF_GPB1.filters.Sk, nan(ny, ny, MKF_GPB1.nh)));
assert(isequaln(MKF_GPB1.filters.Kf, nan(n, ny, MKF_GPB1.nh)));
assert(isequaln(MKF_GPB1.filters.Pk, nan(n, n, MKF_GPB1.nh)));
assert(isequaln(MKF_GPB1.filters.Xk_est, nan(n, 1, MKF_GPB1.nh)));
assert(isequaln(MKF_GPB1.filters.Yk_est, nan(ny, 1, MKF_GPB1.nh)));
assert(isequal(MKF_GPB1.filters.Pkp1, repmat(P0, 1, 1, MKF_GPB1.nh)));
assert(isequal(MKF_GPB1.filters.Xkp1_est, repmat(x0, 1, 1, MKF_GPB1.nh)));

% Also with initial prior probability values
p_seq_g_Yk_init = [0.6; 0.4];
MKF_GPB1x0p0 = MKFObserverGPB1(models,P0,T,"MKF_GPB1x0",x0, ...
    p_seq_g_Yk_init);
assert(isequaln(MKF_GPB1x0p0.x0, x0))
assert(isequaln(MKF_GPB1x0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF_GPB1x0p0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))

% Define GPB2
MKF_GPB2 = MKFObserverGPB2(models,P0,T,"MKF_GPB2",x0);

assert(strcmp(MKF_GPB2.type, "MKF_GPB2"))
assert(isequal(MKF_GPB2.models, models))
assert(isequal(MKF_GPB2.Ts, Ts))
assert(isequal(MKF_GPB2.P0, P0))
assert(isequal(MKF_GPB2.T, T))
assert(strcmp(MKF_GPB2.label, "MKF_GPB2"))
assert(MKF_GPB2.nh == 4)
assert(MKF_GPB2.n == n)
assert(MKF_GPB2.nu == nu)
assert(MKF_GPB2.ny == ny)
assert(MKF_GPB2.nj == 2)
assert(isequal(MKF_GPB2.T, T))
assert(isequal(MKF_GPB2.xkp1_est, x0))
assert(isequal(MKF_GPB2.Pkp1, P0))
assert(isequaln(MKF_GPB2.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB2.Pk, nan(n)))
assert(isequaln(MKF_GPB2.yk_est, nan(ny, 1)))
% assert(isequaln(MKF_GPB2.rkm1, [0 1 0 1]'))
% assert(isequaln(MKF_GPB2.rk, [0 0 1 1]'))
% assert(isequal(MKF_GPB2.rk, r0))
nh = 4;
assert(MKF_GPB2.nh == nh)
assert(isequal(MKF_GPB2.p_seq_g_Yk_init, ones(nh, 1) ./ nh))
assert(isequaln(MKF_GPB2.p_yk_g_seq_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_GPB2.p_rk_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_GPB2.p_seq_g_Ykm1, nan(nh, 1)))
assert(isequaln(MKF_GPB2.filters.Xkp1_est, repmat(x0, 1, 1, nh)))
assert(isequaln(MKF_GPB2.filters.Pkp1, repmat(P0, 1, 1, nh)))
assert(isequaln(MKF_GPB2.filters.Kf, nan(n, ny, nh)))
assert(isequaln(MKF_GPB2.filters.Sk, nan(ny, ny, nh)))
assert(isequaln(MKF_GPB2.merged.Xk_est, nan(n, 1, nj)))
assert(isequaln(MKF_GPB2.merged.Pk, nan(n, n, nj)))
assert(isequaln(MKF_GPB2.merged.Yk_est, nan(ny, 1, nj)))
assert(isequaln(MKF_GPB2.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB2.Pk, nan(n)))
assert(isequaln(MKF_GPB2.yk_est, nan(ny, 1)))

% Run simulation

% Choose observers to include in simulation
%observers = {KF1, KF2, SKF};
%observers = {SKF, MKF_AMM1, MKF_GPB1};
%observers = {MKF_GPB1, MKF_GPB2};
observers = {KF1, KF2, MKF_AMM1, MKF_GPB1, MKF_GPB2, SKF};
%observers = {KF1, KF2, SKF, MKF_AMM1};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = [];

% Simulate observers - without measurement noise (Y)
[Xk_est,Yk_est,DiagPk,MKF_vars] = run_simulation_obs(Y,U,Gamma,seq, ...
    observers,f_mkf);

% Combine and display results
% sim_results1 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs)
% mkf_sim_results = table(t,Gamma,U,X,Y,Ym,MKF_i,MKF_p_seq_g_Yk,MKF_trP_obs)

% % Plot of observer estimates compared to true outputs
% figure(2); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,escape_latex_chars(obs_labels))

% Convert to table
E_obs = array2table(Y - Yk_est, 'VariableNames', obs_labels);
Xk_est = array2table(Xk_est, 'VariableNames', obs_labels);

% Compute mean-squared error
rmses = table(sqrt(mean(E_obs{:,:}.^2))', 'RowNames', obs_labels, ...
    'VariableNames', {'RMSE'});

% Compare to GPB1 & GPB2 simulation results using
% adapted code from Zhao et al.
filename = "GPB_sim_results.csv";
sim_results_GPB = readtable(fullfile(data_dir,filename));

% figure(10); clf;
% 
% ax1 = subplot(3,1,1);
% plot(t, X, t, sim_results_GPB.X)
% xlabel("Time",'Interpreter','latex')
% ylabel("$x(k|k)$",'Interpreter','latex')
% obs_labels = {'MKFObserverGPB1','GPB1_estimation'};
% legend(escape_latex_chars(obs_labels),'Interpreter','latex', ...
%     'location','best')
% set(gca, 'TickLabelInterpreter','latex')
% grid on
% title("True system state",'Interpreter','latex')
% 
% ax2 = subplot(3,1,2);
% plot(t, Xk_est{:,'MKF_GPB1'}, t, sim_results_GPB.Xk_est_GPB1)
% xlabel("Time",'Interpreter','latex')
% ylabel("$\hat{x}(k|k)$",'Interpreter','latex')
% obs_labels = {'MKFObserverGPB1','GPB1_estimation'};
% legend(escape_latex_chars(obs_labels),'Interpreter','latex', ...
%     'location','best')
% set(gca, 'TickLabelInterpreter','latex')
% grid on
% title("Observer estimates",'Interpreter','latex')
% 
% ax3 = subplot(3,1,3);
% plot(t, Xk_est{:,'MKF_GPB2'}, t, sim_results_GPB.Xk_est_GPB2)
% xlabel("Time",'Interpreter','latex')
% ylabel("$\hat{x}(k|k)$",'Interpreter','latex')
% ylim([-1 11])
% obs_labels = {'MKFObserverGPB2','GPB2_estimation'};
% legend(escape_latex_chars(obs_labels),'Interpreter','latex', ...
%     'location','best')
% set(gca, 'TickLabelInterpreter','latex')
% grid on
% title("Observer estimates",'Interpreter','latex')
% 
% linkaxes([ax1 ax2 ax3],'xy')

% Note: Set measurement noise to 0.1 for these tests.

% Check KF1 was accurate before system switched
assert(mean(abs(E_obs{t <= 9.5, "KF1"})) < 0.1)  % This was 1e-5 before is that possible?

% Check KF2 was accurate after system switched
assert(max(abs(E_obs{t > 20, "KF2"})) < 0.1)

if ismember("SKF_S",obs_labels)
    % Check SKF matches KF1 before system switched
    assert(max(abs(Xk_est{t <= 9.5, "KF1"} - Xk_est{t <= 9.5, "SKF_S"})) < 1e-4)
    % Check SKF converges to KF2 after system switched
    assert(max(abs(Xk_est{t > 20, "KF2"} - Xk_est{t > 20, "SKF_S"})) < 0.01)
else
    warning("SKF_S not included in simulation")
end

if ismember("MKF_AMM1",obs_labels)
    % Check MKF_AMM converged to KF1 before system switched - it does not
    %assert(max(abs(Xk_est{t > 5 & t <= 9.5, "KF1"} - Xk_est{t > 5 & t <= 9.5, "MKF_AMM1"})) < 1e-10)
    % Check MKF_AMM converged to KF2 after system switched - it does not
    % TODO: Should it switch?
    %assert(max(abs(Xk_est{t > 20, "KF2"} - Xk_est{t > 20, "MKF_AMM1"})) < 1e-10)
else
    warning("MKF_AMM1 not included in simulation")
end

% Check GPB1
if any(strcmp("MKF_GPB1", Xk_est.Properties.VariableNames))
    % Check GPB1 is close to KF1 before system switched
    assert(abs(Xk_est{t == 9.5, "KF1"} - Xk_est{t == 9.5, "MKF_GPB1"}) < 1e-3)
    % Check GPB1 is close to KF2 after system switched
    assert(max(abs(Xk_est{t > 20, "KF2"} - Xk_est{t > 20, "MKF_GPB1"})) < 0.005)
end

% Check root-mean-squared errors
rmse_test = table([1.9509 0.519 1.9509 0.0310 0.0310 0.0000]', ...
    'RowNames', ["KF1"  "KF2" "MKF_AMM1" "MKF_GPB1" "MKF_GPB2" "SKF_S"], ...
    'VariableNames', {'RMSE'} ...
);
for i = 1:size(rmses, 1)
    label = rmses.Properties.RowNames{i};
    assert(round(rmses{label, 'RMSE'}, 4) == rmse_test{label, 'RMSE'})
end

% Test reset methods
% Should set observer variables to original initial conditions
KF1.reset()
KF2.reset()
MKF_AMM1.reset()
MKF_GPB1.reset()
SKF.reset()

assert(isequal(MKF_AMM1.xkp1_est, MKF_AMM1.x0))
assert(isequal(MKF_AMM1.Pkp1, MKF_AMM1.P0))
assert(isequal(MKF_AMM1.p_seq_g_Yk, MKF_AMM1.p_seq_g_Yk_init))

assert(isequal(MKF_GPB1.xkp1_est, MKF_GPB1.x0))
assert(isequal(MKF_GPB1.Pkp1, MKF_GPB1.P0))
assert(isequal(MKF_GPB1.p_seq_g_Yk, MKF_GPB1.p_seq_g_Yk_init))
assert(isequaln(MKF_GPB1.p_yk_g_seq_Ykm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.p_rk_g_Ykm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.p_rk_g_rkm1, nan(MKF_GPB1.nh, 1)))  % TODO: Is this correct?
assert(isequaln(MKF_GPB1.p_seq_g_Ykm1, nan(MKF_GPB1.nh, 1)))

% Redefine a new observer (identical to above)
MKF_GPB1_new = MKFObserverGPB1(models,P0,T,'MKF_GPB1',x0);
assert(isequaln(MKF_GPB1_new, MKF_GPB1))
MKF_GPB1_new.label = "MKF_GPB1_new";

% Make a copy
MKF_GPB1_copy = MKF_GPB1_new.copy();
assert(isequaln(MKF_GPB1_copy, MKF_GPB1_new))
MKF_GPB1_copy.label = "MKF_GPB1_copy";

% Choose observers to include in simulation
observers = {KF1, KF2, SKF, MKF_GPB1, MKF_GPB1_new, MKF_GPB1_copy, MKF_GPB2};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = [];

% Simulate observers - with measurement noise (Ym)
[Xk_est,Yk_est,DiagPk,MKF_vars] = run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf);

% % Combine and display results
% sim_results = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs);
% writetable(sim_results, "results/test_MKFO_sim_results.csv");

% % Display results from MKF observer
% sim_results_MKF = [ ...
%     table(t) ... 
%     table(MKF_vars.K_obs) ...
%     table(MKF_vars.trP_obs) ...
%     table(MKF_vars.i) ...
%     table(MKF_vars.p_seq_g_Yk) ...
% ];
% writetable(sim_results_MKF, "results/test_MKFO_sim_results_MKF.csv");

% Check final state estimates
test_X_est = [-1.2111    9.9249    9.9249    9.9249    9.9249    9.9249    9.9249];
assert(isequal(round(Xk_est(t == t(end), :), 4), test_X_est))

% Check final error covariance estimates
test_DiagP = [    0.0150    0.0225    0.0225    0.0225    0.0225    0.0225    0.0225];
assert(isequal(round(DiagPk(t == t(end), :), 4), test_DiagP))

% % Plot of observer estimates compared to true outputs
% figure(3); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,escape_latex_chars(obs_labels))

% Convert to table
E_obs = array2table(Y - Yk_est, 'VariableNames', obs_labels);
Xk_est = array2table(Xk_est, 'VariableNames', obs_labels);

% Compute mean-squared error
rmses = table(sqrt(mean(E_obs{:,:}.^2))', 'RowNames', obs_labels, ...
    'VariableNames', {'RMSE'});

% Check root-mean-squared errors
rmse_test = table([1.9509 0.496 0.0302 0.0471 0.0471 0.0471 0.0448]', ...
    'RowNames', ["KF1" "KF2" "SKF_S" "MKF_GPB1" "MKF_GPB1_new" "MKF_GPB1_copy" "MKF_GPB2"], ...
    'VariableNames', {'RMSE'} ...
);
for i = 1:size(rmses, 1)
    label = rmses.Properties.RowNames{i};
    assert(round(rmses{label, 'RMSE'}, 4) == rmse_test{label, 'RMSE'})
end

% % Plot selected observers
% figure(4); clf
% plot_obs_estimates(t,X,X_est(:,[3 4]),Y,Y_est(:,[3 4]),obs_labels([3 4]))


%% Simulation test on 2x2 system with random shock

% Sample time
Ts = 1;

% NOTE: this is a previous version of the system with lower
% coupling (-0.2) and epsilon = [0.01; 0.01].

% Discrete time state space model
A = [ 0.8890       0     1 -0.2;
           0  0.8890  -0.2    1;
           0       0     1    0;
           0       0     0    1];
B = [    1 -0.2  0  0;
      -0.2    1  0  0;
         0    0  1  0;
         0    0  0  1];
C = [ 0.1110       0  0  0;
           0  0.1110  0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Designate measured input and output signals
u_meas = [true; true; false; false];
y_meas = [true; true];

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.1; 0.1];
sigma_wp = [0.01 1; 0.01 1];

% Simulation settings
nT = 100;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = [5 10];
du0 = [1; 1];
% When you make the shock larger the MKF observers
% do better
%du0 = [2; 2];

% Measured input
%U = (idinput(size(t)) + 1)/2;
U = zeros(nT+1, 2);
U(t >= 1, 1) = -1;

% Disturbance input
% This is used by the SKF observer
Gamma = zeros(nT+1, 2);
Gamma(t == t_shock(1), 1) = 1;
Gamma(t == t_shock(2), 2) = 1;
Wp = du0' .* Gamma;

U_sim = [U Wp];

% Observer parameters (same for all observers)
P0 = 10000*eye(n);
x0 = [0.5 -0.1 0.1 -0.1]';  % optional
y0 = C * x0;  % optional

% 4 possible process noise covariance matrices
Qj = {diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]), ...
      diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]), ...
      diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2]), ...
      diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,2)^2])};
R = diag(sigma_M.^2);

% 3 system models
m1 = struct();
m1.A = A;
m1.B = Bu;
m1.C = C;
m1.Q = Qj{1};
m1.R = R;
m1.Ts = Ts;
m2 = struct();
m2.A = A;
m2.B = Bu;
m2.C = C;
m2.Q = Qj{2};
m2.R = R;
m2.Ts = Ts;
m3 = struct();
m3.A = A;
m3.B = Bu;
m3.C = C;
m3.Q = Qj{3};
m3.R = R;
m3.Ts = Ts;
models = {m1, m2, m3};
[nj, n, nu, ny, Ts] = check_models(models);

% Custom MKF observer with correct shock sequence
seq = repmat({ones(1, nT+1)}, 4, 1);
seq{2}(t == t_shock(1)) = 2;  % shock 1
seq{3}(t == t_shock(2)) = 3;  % shock 2
seq{4}(t == t_shock(1)) = 2;  % both
seq{4}(t == t_shock(2)) = 3;
p_rk = [1-epsilon epsilon]';
Z = [1 1; 2 1; 1 2];  % combinations
p_rk = prod(prob_rk(Z', p_rk), 1)';
p_rk = p_rk ./ sum(p_rk);  % normalized
T = repmat(p_rk', nj, 1);

% Define MKF observers

% 1. MKF_AMM - makes no sense for RODD disturbance
MKF_AMM1 = MKFObserverAMM(models,P0,'MKF_AMM1');

assert(strcmp(MKF_AMM1.type, "MKF_AMM"))
assert(isequal(MKF_AMM1.models, models))
assert(isequal(MKF_AMM1.Ts, Ts))
assert(isequal(MKF_AMM1.P0, P0))
assert(strcmp(MKF_AMM1.label, "MKF_AMM1"))
assert(MKF_AMM1.nh == 3)
assert(MKF_AMM1.n == n)
assert(MKF_AMM1.nu == nu)
assert(MKF_AMM1.ny == ny)
assert(MKF_AMM1.nj == 3)
assert(isequal(MKF_AMM1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_AMM1.Pkp1, P0))
assert(isequaln(MKF_AMM1.xk_est, nan(n, 1)))
assert(isequaln(MKF_AMM1.Pk, nan(n)))
assert(isequaln(MKF_AMM1.yk_est, nan(ny, 1)))
assert(isequaln(MKF_AMM1.x0, zeros(n, 1)))
assert(isequal(MKF_AMM1.p_seq_g_Yk_init, ones(3, 1) ./ 3))

% Redefine this time with initial conditions
MKF_AMM1x0 = MKFObserverAMM(models,P0,'MKF_AMM1x0',x0);
assert(isequaln(MKF_AMM1x0.x0, x0))
assert(isequal(MKF_AMM1x0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0.Pkp1, P0))

% Also with initial prior probability values
p_seq_g_Yk_init = [0.6; 0.4];
MKF_AMM1x0p0 = MKFObserverAMM(models,P0,"MKF_GPB1x0",x0, ...
    p_seq_g_Yk_init);
assert(isequaln(MKF_AMM1x0p0.x0, x0))
assert(isequaln(MKF_AMM1x0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF_AMM1x0p0.xkp1_est, x0))
assert(isequal(MKF_AMM1x0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))

% 2. MKF_GPB1
MKF_GPB1 = MKFObserverGPB1(models,P0,T,'MKF_GPB1');

assert(strcmp(MKF_GPB1.type, "MKF_GPB1"))
assert(isequal(MKF_GPB1.models, models))
assert(isequal(MKF_GPB1.Ts, Ts))
assert(isequal(MKF_GPB1.P0, P0))
assert(isequal(MKF_GPB1.T, T))
assert(strcmp(MKF_GPB1.label, "MKF_GPB1"))
assert(MKF_GPB1.nh == 3)
assert(MKF_GPB1.n == n)
assert(MKF_GPB1.nu == nu)
assert(MKF_GPB1.ny == ny)
assert(MKF_GPB1.nj == 3)
assert(isequal(MKF_GPB1.T, T))
assert(isequal(MKF_GPB1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_GPB1.Pkp1, P0))
assert(isequaln(MKF_GPB1.xk_est, nan(n, 1)))
assert(isequaln(MKF_GPB1.Pk, nan(n)))
assert(isequaln(MKF_GPB1.yk_est, nan(ny, 1)))
assert(isequaln(MKF_GPB1.p_yk_g_seq_Ykm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.p_rk_g_Ykm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.p_rk_g_rkm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.p_seq_g_Ykm1, nan(MKF_GPB1.nh, 1)))
assert(isequaln(MKF_GPB1.x0, zeros(n, 1)))
assert(isequal(MKF_GPB1.p_seq_g_Yk_init, ones(3, 1) ./ 3))
assert(isequaln(MKF_GPB1.filters.Sk, nan(ny, ny, MKF_GPB1.nh)));
assert(isequaln(MKF_GPB1.filters.Kf, nan(n, ny, MKF_GPB1.nh)));
assert(isequaln(MKF_GPB1.filters.Pk, nan(n, n, MKF_GPB1.nh)));
assert(isequaln(MKF_GPB1.filters.Xk_est, nan(n, 1, MKF_GPB1.nh)));
assert(isequaln(MKF_GPB1.filters.Yk_est, nan(ny, 1, MKF_GPB1.nh)));
assert(isequal(MKF_GPB1.filters.Pkp1, repmat(P0, 1, 1, MKF_GPB1.nh)));
assert(isequal(MKF_GPB1.filters.Xkp1_est, zeros(n, 1, MKF_GPB1.nh)));

% Redefine this time with initial conditions
MKF_GPB1x0 = MKFObserverGPB1(models,P0,T,"MKF_GPB1x0",x0);
assert(isequaln(MKF_GPB1x0.x0, x0))
assert(isequal(MKF_GPB1x0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0.Pkp1, P0))
assert(isequal(MKF_GPB1x0.p_seq_g_Yk_init, ones(3, 1) ./ 3))

% Also with initial prior probability values
p_seq_g_Yk_init = [0.6; 0.4];
MKF_GPB1x0p0 = MKFObserverGPB1(models,P0,T,"MKF_GPB1x0",x0, ...
    p_seq_g_Yk_init);
assert(isequaln(MKF_GPB1x0p0.x0, x0))
assert(isequaln(MKF_GPB1x0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))
assert(isequal(MKF_GPB1x0p0.xkp1_est, x0))
assert(isequal(MKF_GPB1x0p0.p_seq_g_Yk_init, p_seq_g_Yk_init))

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
Aj = repmat({A}, 1, nj);
Buj = repmat({Bu}, 1, nj);
Cj = repmat({C}, 1, nj);
Rj = repmat({R}, 1, nj);
SKF = SKFObserverS(models,P0,seq{4},"SKF",x0);

% Choose observers to test
observers = {MKF_AMM1, MKF_GPB1, SKF};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Simulate system
X = zeros(nT+1,n);
Y = zeros(nT+1,ny);
xk = zeros(n,1);

for i = 1:nT+1

    % Inputs
    uk = U_sim(i,:)';

    % Compute y(k)
    yk = C*xk + D*uk;

    % Store results
    X(i, :) = xk';
    Y(i, :) = yk';

    % Compute x(k+1)
    xk = A*xk + B*uk;

end

% Check simulation output is correct
[Y2, t, X2] = lsim(Gpss, U_sim, t);
assert(isequal(X, X2))
assert(isequal(Y, Y2))

% Choose measurement noise for plant
sigma_MP = [0; 0];  % Set to zero for testing
Ym = Y + sigma_MP'.*randn(nT+1, ny);

% Identify which observer to log MKF data for
f_mkf = 1;

% Simulate observers
[Xk_est,Yk_est,DiagPk,MKF_vars] = run_simulation_obs(Ym,U,Gamma,seq,observers,f_mkf);

% Output estimation errors
E_obs = repmat(Y, 1, n_obs) - Yk_est;

% Combine and display results
sim_results1 = table(t,Gamma,U,X,Y,Ym,Xk_est,Yk_est,E_obs);

% Plot observer estimates
% figure(5); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,obs_labels)

% Average results and convert to table
E_obs_avg = array2table(reshape(mean(reshape(E_obs,[],ny,n_obs), 2),[],n_obs), ...
    'VariableNames', obs_labels);
Xk_est_avg = array2table(reshape(mean(reshape(Xk_est,[],n,n_obs), 2),[],n_obs), ...
    'VariableNames', obs_labels);
DiagPk_avg = array2table(reshape(mean(reshape(DiagPk,[],n,n_obs), 2),[],n_obs), ...
    'VariableNames', obs_labels);

% Check final state estimates
test_X_est_avg = [2.293238    2.301738    2.301600];
assert(isequal(round(Xk_est_avg{t == t(end), :}, 6), test_X_est_avg))

% Check final error covariance estimates
% TODO: Haven't checked if these are correct.
test_DiagPk_avg = [0.0427    0.1547    0.0427];
assert(isequal(round(DiagPk_avg{t == t(end), :}, 4), test_DiagPk_avg))

% Compute mean-squared error
rmses = array2table(reshape(sqrt(mean(reshape(E_obs,[],ny,n_obs).^2, 2)),[],n_obs), ...
    'VariableNames', obs_labels);
rmses = array2table(round(mean(rmses{:,:}, 1), 6), 'VariableNames', obs_labels);
%fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

% % Display results of last simulation
% 
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% K_obs = sim_results.K_obs;
% trP_obs = sim_results.trP_obs;
% 
% table(t,alpha,U,Wp,X,Y,Y_m,X_est,Y_est,E_obs)
% 
% % Display gains and trace of covariance matrix
% table(t, cell2mat(K_obs), cell2mat(trP_obs), ...
%     'VariableNames', {'t', 'K{1}, K{2}', 'trace(P{1}), trace(P{2})'})
% 
% % Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'})

% Results
MSE_test_values = array2table(...
    [0.025273 0.010339 0.002565], ...
    'VariableNames', {'MKF_AMM1', 'MKF_GPB1', 'SKF'});
assert(isequal(rmses, MSE_test_values))


%% Test copy methods

% Define system

% Sample period
Ts = 0.5;

% Discrete time state space models
% Model #1
model1 = struct();
model1.A = 0.7;
model1.B = 1;
model1.C = 0.3;
model1.D = 0;
model1.Ts = Ts;
%Gpss1 = ss(A1,B1,C1,D1,Ts);

% Model #2
model2 = struct();
model2.A = 0.9;
model2.B = 1;
model2.C = -0.3;  % -ve gain!
model2.D = 0;
model2.Ts = Ts;
%Gpss2 = ss(A2,B2,C2,D2,Ts);

% Switching system
models = {model1, model2};

% Check dimensions
[nj, n, nu, ny, Ts] = check_models(models);

% Observer parameters (same for all observers)
P0 = 10000;
x0 = 0.5;  % optional
models{1}.Q = 0.01;
models{1}.R = 0.1^2;
models{2}.Q = 0.01;
models{2}.R = 0.1^2;

% Transition probabilities
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
assert(all(sum(T, 2) == 1))

% System indicator sequences
nT = 100;
seq1 = {
    zeros(1, nT+1);
    [zeros(1, 20) ones(1, nT+1-20)];  % equal to Gamma'
    [zeros(1, 40) ones(1, nT+1-40)];
    ones(1, nT+1);
 };

% Define MKF observers

% Define multi-model observer with initial conditions
MKF = MKFObserverGPB1(models,P0,T,'MKF',x0);

% Test handle copy
MKF_hcopy = MKF;
assert(isequaln(MKF_hcopy, MKF))  % same values
assert(MKF_hcopy == MKF)  % must be same object

MKF.x0 = 1.0;
assert(isequal(MKF_hcopy.x0, 1.0))

% Test true copy
MKF_copy = MKF.copy();
assert(isequaln(MKF_copy, MKF))  % same values
assert(MKF_copy ~= MKF)  % must not be same object

MKF.label = "New name";
assert(~isequal(MKF_copy.label, "New name"))


% THIS HAS BEEN REPLACED BY run_simulation_obs.m
%
% function [Xk_est,Yk_est,DiagPk,Xkp1_est,Ykp1_est,DiagPkp1,MKF_K_obs,MKF_trP_obs, ...
%     MKF_i,MKF_p_seq_g_Yk] = run_simulation_obs(Ym,U,observers,f_mkf)
% % Simulate observers
% 
%     nT = size(Ym, 1) - 1;
%     ny = size(Ym, 2);
%     n_obs = numel(observers);
%     n = size(observers{1}.xkp1_est, 1);
% 
%     obs_mkf = observers{f_mkf};
%     n_filters = obs_mkf.nh;
% 
%     Xk_est = zeros(nT+1, n*n_obs);
%     Yk_est = zeros(nT+1, ny*n_obs);
%     Xkp1_est = zeros(nT+1, n*n_obs);
%     Ykp1_est = zeros(nT+1, ny*n_obs);
%     DiagPk = nan(nT+1, n*n_obs);
%     DiagPkp1 = nan(nT+1, n*n_obs);
%     MKF_K_obs = cell(nT+1, n*n_filters);
%     MKF_trP_obs = nan(nT+1, n_filters);
%     MKF_i = nan(nT+1, 2);
%     MKF_p_seq_g_Yk = nan(nT+1, n_filters);
% 
%     for i = 1:nT+1
% 
%         yk = Ym(i, :)';
%         uk = U(i, :)';
% 
%         % Arrays to store results
%         xk_est = nan(1, n*n_obs);
%         yk_est = nan(1, ny*n_obs);
%         xkp1_est = nan(1, n*n_obs);
%         ykp1_est = nan(1, ny*n_obs);
%         diagPk = nan(1, n*n_obs);
%         diagPkp1 = nan(1, n*n_obs);
% 
%         % Update observers
%         for f = 1:n_obs
%             obs = observers{f};
%             obs.update(yk, uk);
%             if f == f_mkf
%                 if isprop(obs, "i")
%                     MKF_i(i, :) = obs.i;
%                 else
%                     MKF_i(i, :) = nan;
%                 end
%                 MKF_p_seq_g_Yk(i, :) = obs.p_seq_g_Yk';
%                 for j = 1:obs.nh
%                     switch obs.type
%                         case {"KF", "SKF", "MKF"}
%                             % Only works for ny = 1
%                             %MKF_K_obs{i, j} = obs.filters{j}.K';
%                             MKF_trP_obs(i, j) = trace(obs.filters{j}.P);
%                         case {"KFF", "SKFF"}
%                             %MKF_K_obs{i, j} = obs.filters{j}.Kf';
%                             MKF_trP_obs(i, j) = trace(obs.filters{j}.Pkp1);
%                         case {"MKF_AMM", "MKF_GPB1"}
%                             MKF_trP_obs(i, j) = trace(obs.Pkf(:, :, j));
%                     end
%                 end
%                 switch obs.type
%                 end
%             end
%             if isprop(obs,'xk_est')
%                 xk_est(1, (f-1)*n+1:f*n) = obs.xk_est';
%                 yk_est(1, (f-1)*ny+1:f*ny) = obs.yk_est';
%                 diagPk(1, (f-1)*n+1:f*n) = diag(obs.Pk)';
%             end
%             if isprop(obs,'xkp1_est')
%                 xkp1_est(1, (f-1)*n+1:f*n) = obs.xkp1_est';
%                 ykp1_est(1, (f-1)*ny+1:f*ny) = obs.ykp1_est';
%                 diagPkp1(1, (f-1)*n+1:f*n) = diag(obs.Pkp1)';
%             end
%         end
% 
%         % Record observer estimates
%         Xk_est(i, :) = xk_est;
%         Yk_est(i, :) = yk_est;
%         Xkp1_est(i, :) = xkp1_est;
%         Ykp1_est(i, :) = ykp1_est;
%         DiagPk(i, :) = diagPk;
%         DiagPkp1(i, :) = diagPkp1;
% 
%     end
% end
