% Test KalmanFilter, KalmanFilterF, and KalmanFilterSS classes
% (See class definitions in KalmanFilter.m, KalmanFilterF.m,
% and KalmanFilterSS.m)

clear all

addpath("../plot-utils/")

% Folder containing test data
results_dir = 'results';


%% SISO system example from GEL-7029
% See file Kalman_Filter.mlx

sys_test_siso

% Initial state of system
x0 = [0.1; 0.5];

% Check if benchmark simulation data file exists
filename = 'KF_sim_benchmark.csv';
if ~isfile(fullfile(results_dir, filename))
    error("Test data file '%s' not found.\n Run 'Kalman_Filter_benchmark.mlx' to generate benchmark data.", filename)
end
bench_sim_results = readtable(fullfile(results_dir, filename));

% Kalman filter parameters (not same as those of system)
Q = diag([0.1; 0.1]);
R = 0.5;
N = zeros(n,ny);
x0_est = zeros(2, 1);
% Note: Qp and Rp are system noise variances

% Prepare a struct of model parameters
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.Ts = Ts;
model.Q = Q;
model.R = R;

% Steady-state Kalman filter in prediction form
label = "KFPSS";
KFPSS = KalmanFilterPSS(model,label,x0_est);
assert(strcmp(KFPSS.type, "KFPSS"))
assert(isequal(KFPSS.model, model))
K_calc = A * KFPSS.Pkp1 * C' * (C * KFPSS.Pkp1 * C' + R)^-1;
assert(max(abs(KFPSS.K - K_calc)) < 1e-12)
assert(isequal(round(KFPSS.K, 6), [0.772750; 0.755731]))
assert(isequal(round(KFPSS.Pkp1, 6), [1.509786 1.216953; 1.216953 1.219071]))
assert(isequal(KFPSS.label, label))
assert(isequal(KFPSS.xkp1_est, x0_est))
assert(KFPSS.ykp1_est == C*x0_est)
assert(KFPSS.n == n)
assert(KFPSS.nu == nu)
assert(KFPSS.ny == ny)

% Steady-state Kalman filter in filtering form
label = "KFFSS";
KFFSS = KalmanFilterFSS(model,label,x0_est);
assert(strcmp(KFFSS.type, "KFFSS"))
assert(isequal(KFFSS.model, model))
Kf_calc = KFFSS.Pkp1 * C' * (C * KFFSS.Pkp1 * C' + R)^-1;
assert(max(abs(KFFSS.Kf - Kf_calc)) < 1e-12)
assert(isequal(round(KFFSS.Kf, 6), [0.755731; 0.669133]))
assert(isequal(round(KFFSS.Pkp1, 6), [1.509786 1.216953; 1.216953 1.219071]))
assert(isequal(KFFSS.label, label))
assert(isequal(KFFSS.xkp1_est, x0_est))
assert(KFFSS.ykp1_est == C*x0_est)
assert(isequaln(KFFSS.xk_est, nan(n, 1)))
assert(isequaln(KFFSS.yk_est, nan))
assert(KFFSS.n == n)
assert(KFFSS.nu == nu)
assert(KFFSS.ny == ny)

% Re-define with no initial state specified (should be set to zero)
KFPSSx0 = KalmanFilterPSS(model,"KFPSSx0");
assert(isequal(round(KFPSS.K, 6), round(KFPSSx0.K, 6)))
assert(isequal(round(KFPSS.Pkp1, 6), round(KFPSSx0.Pkp1, 6)))
assert(isequal(KFPSSx0.xkp1_est, zeros(n, 1)))
assert(KFPSSx0.ykp1_est == 0)
KFFSSx0 = KalmanFilterFSS(model,"KFFSSx0");
assert(isequal(round(KFFSS.Kf, 6), round(KFFSSx0.Kf, 6)))
assert(isequal(round(KFFSS.Pkp1, 6), round(KFFSSx0.Pkp1, 6)))
assert(isequal(KFFSSx0.xkp1_est, zeros(n, 1)))
assert(KFFSSx0.ykp1_est == 0)

% Define dynamic Kalman filters with KalmanFilterP
% and KalmanFilterF classes
P0 = diag([1e-4 1e-4]);

KF_old = KalmanFilter(A,B,C,Ts,P0,Q,R,"KF_old",x0_est);
assert(isequal(KF_old.A, A))
assert(isequal(KF_old.B, B))
assert(isequal(KF_old.C, C))
assert(isequal(KF_old.Ts, Ts))
assert(isequal(KF_old.Pkp1, P0))
assert(isequal(KF_old.Q, Q))
assert(isequal(KF_old.R, R))
assert(all(isnan(KF_old.K)))
assert(isequal(KF_old.xkp1_est, x0_est))
assert(KF_old.ykp1_est == KF_old.C * x0_est)

% New version in filtering form
KFF = KalmanFilterF(model,P0,"KFF",x0_est);
assert(isequal(KFF.model, model))
assert(all(isnan(KFF.Kf)))
assert(isequal(KFF.xkp1_est, x0_est))
assert(KFF.ykp1_est == KFF.model.C * x0_est)
assert(isequal(KFF.Pkp1, P0))
assert(all(isnan(KFF.xk_est)))
assert(all(isnan(KFF.yk_est)))
assert(all(isequaln(KFF.Pk, nan(2))))

% New version in prediction form
KFP = KalmanFilterP(model,P0,"KFP",x0_est);
assert(strcmp(KFP.type, "KFP"))
assert(isequal(KFP.model, model))
assert(all(isnan(KFP.K)))
assert(isequal(KFP.xkp1_est, x0_est))
assert(KFP.ykp1_est == KFP.model.C * x0_est)
assert(isequal(KFP.Pkp1, P0))
assert(all(isnan(KFP.xk_est)))
assert(all(isnan(KFP.yk_est)))

% number of points to simulate
nT = 100;

% NOTE: To match benchmark results, rng seed and initialization
% of input data must match code in Kalman_Filter_benchmark.mlx
% rng(0)  % seed random number generator
% v=sqrt(Rp)*randn(Np,1); % measurement noise for the whole simulation
% w=sqrt(Qp)*randn(2,Np); % proces noise for the whole simulation

% seed random number generator
rng(0)

% Measurement noise for the whole simulation
V = sqrt(Rp)*randn(1, nT)';

% Process noise for the whole simulation
W = (sqrt(Qp)*randn(2, nT))';

u0 = 1;  % initial value of u
x0 = (eye(length(A)) - A) \ (B*u0);  % steady-state value of x

% Intialize system (at k = 0)
xk = x0;

% Input signal
U = [zeros(10,1); ones(nT-10, 1)];  % process input for the whole simulation

% Observers to simulate
observers = {KFPSS, KFFSS, KF_old, KFP, KFF};
obs_labels = cellfun(@(obs) obs.label, observers);
n_obs = numel(observers);

% Matrices to collect simulation data
Xk = zeros(nT, n); % process states
Yk = zeros(nT, ny); % process states
Xk_est = cell(1, n_obs);  % process state estimates
Yk_est = cell(1, n_obs);  % process output estimates
Xkp1_est = cell(1, n_obs);
Ykp1_est = cell(1, n_obs);
for f = 1:n_obs
    Xk_est{f} = nan(nT, n);
    Yk_est{f} = nan(nT, ny);
    Xkp1_est{f} = nan(nT, n);
    Ykp1_est{f} = nan(nT, ny);
end

% Time vector
t = Ts * (0:nT-1)';
for i = 1:nT

    % Inputs
    vk = V(i, :)';
    uk = U(i, :)';
    wk = W(i, :)';

    % Process output and input in current timestep
    yk = C*xk + vk;

    % Record process states and output
    Xk(i, :) = xk';
    Yk(i, :) = yk';

    % Check predictions in next time step
    % KFF.predict();
    % assert(all(abs(A * KFF.xk_est + B * U(i) - KFF.xkp1_est) < 1e-14))

    for f = 1:n_obs

        obs = observers{f};

        % Record Kalman filter predictions
        if isprop(obs, "xkp1_est") || isfield(obs, "xkp1_est")
            Xkp1_est{f}(i, :) = obs.xkp1_est';
        end
        if isprop(obs, "ykp1_est") || isfield(obs, "ykp1_est")
            Ykp1_est{f}(i, :) = obs.ykp1_est';
            assert(isequaln(C * obs.xkp1_est, obs.ykp1_est))
        end

        % Update Kalman filter
        if strcmp(obs.type, "KFSS")
            obs = update_KF(obs, uk, yk);
        else
            obs.update(yk, uk);
        end

        % Record Kalman filter updated estimates
        if isprop(obs, "xk_est") || isfield(obs, "xk_est")
            Xk_est{f}(i, :) = obs.xk_est';
        end
        if isprop(obs, "yk_est") || isfield(obs, "yk_est")
            Yk_est{f}(i, :) = obs.yk_est';
            assert(isequaln(C * obs.xk_est, obs.yk_est))
        end

        % TODO: Remove this once not using structs any more
        observers{f} = obs;

    end

    % Process states in next timestep
    xk = A*xk + B*uk + wk;

    assert(all(abs(observers{1}.xkp1_est - KFPSS.xkp1_est) < 1e-12))
    assert(all(abs(observers{1}.ykp1_est - KFPSS.ykp1_est) < 1e-12))

    % Check old and new predictions are the same
    assert(all(abs(KF_old.xkp1_est - KFP.xkp1_est) < 1e-12))
    assert(all(abs(KF_old.ykp1_est - KFP.ykp1_est) < 1e-12))
    assert(all(abs(KF_old.xkp1_est - KFF.xkp1_est) < 1e-12))
    assert(all(abs(KF_old.ykp1_est - KFF.ykp1_est) < 1e-12))

end

% Check system simulation is correct
% Add process noise inputs to system model
Gdss_sim = Gdss;
Gdss_sim.B = [Gdss_sim.B eye(2)];
U_sim = [U(1:100, :) W(1:100, :)];
t_sim = t(1:100);
[Y_sim,t_sim,X_sim] = lsim(Gdss_sim,U_sim,t_sim,x0);
Y_sim = Y_sim + V(1:100);
assert(max(abs(Yk - Y_sim)) < 1e-14)

% % Plot results
%
% figure(1); clf
% labels = ["Process" escape_latex_chars(obs_labels)];
%
% ax1 = subplot(411);
% plot(t, Yk, 'k'); hold on
% for f = 1:n_obs
%     plot(t, Yk_est{f}, 'Linewidth', 2)
% end
% legend(labels, 'Interpreter', 'Latex')
% ylabel('$y_1$', 'Interpreter', 'Latex')
% grid on
% 
% ax2 = subplot(412);
% plot(t, Xk(:, 1), 'k'); hold on
% for f = 1:n_obs
%     plot(t, Xk_est{f}(:, 1), 'Linewidth', 2)
% end
% legend(labels, 'Interpreter', 'Latex')
% ylabel('$x_1$', 'Interpreter', 'Latex')
% grid on
% 
% ax3 = subplot(413);
% plot(t, Xk(:, 2), 'k'); hold on
% for f = 1:n_obs
%     plot(t, Xk_est{f}(:, 2), 'Linewidth', 2)
% end
% legend(labels, 'Interpreter', 'Latex')
% ylabel('$x_2$', 'Interpreter', 'Latex')
% grid on
% 
% ax4 = subplot(414);
% stairs(t', U', 'Linewidth', 2);
% xlabel('Time [s]', 'Interpreter', 'Latex');
% ylabel('$u_1$', 'Interpreter', 'Latex')
% grid on
% 
% linkaxes([ax1 ax2 ax3 ax4], 'x')

% Display results
sim_results = [
    table(t,U) ...
    array2table(Xk, 'VariableNames', {'x1', 'x2'}) ...
    array2table(Yk, 'VariableNames', {'y'})
];
for f = 1:n_obs
    obs = observers{f};
    labels = compose("xk_est_%d_", 1:n) + obs.label;
    sim_results(:, labels) = array2table(Xk_est{f});
end
for f = 1:n_obs
    obs = observers{f};
    labels = compose("xkp1_est_%d_", 1:n) + obs.label;
    sim_results(:, labels) = array2table(Xkp1_est{f});
end
for f = 1:n_obs
    obs = observers{f};
    labels = compose("yk_est_%d_", 1:ny) + obs.label;
    sim_results(:, labels) = array2table(Yk_est{f});
end
for f = 1:n_obs
    obs = observers{f};
    labels = compose("ykp1_est_%d_", 1:ny) + obs.label;
    sim_results(:, labels) = array2table(Ykp1_est{f});
end
%head(sim_results)

% Verify results by comparing with outputs of Kalman_Filter.mlx
%head(bench_sim_results)

% System states
assert(isequal( ...
    round(sim_results{1:100, {'x1', 'x2'}}, 7), ...
    round(bench_sim_results{1:100, {'xNprocess_1', 'xNprocess_2'}}, 7) ...
))

% Check steady-state correction gains compared to that from
% benchmark simulation. See script Kalman_Filter_benchmark.mlx
% >> fprintf("%.6f %.6f\n", KalmanFilter.B')
% 0.250000 0.772750
% 0.000000 0.755731
%
assert(isequal(round(KFPSS.K, 6), [0.772750; 0.755731]))
assert(isequal(round(A * KFFSS.Kf, 6), [0.772750; 0.755731]))

% Predictions of steady-state estimators
assert(isequal( ...
    round(sim_results{:, {'xkp1_est_1_KFPSS', 'xkp1_est_2_KFPSS'}}, 7), ...
    round(bench_sim_results{:, {'xNkalman_1', 'xNkalman_2'}}, 7) ...
))
assert(isequal( ...
    round(sim_results{:, {'xkp1_est_1_KFFSS', 'xkp1_est_2_KFFSS'}}, 7), ...
    round(bench_sim_results{:, {'xNkalman_1', 'xNkalman_2'}}, 7) ...
))

% Calculate state estimation errors and mean-squared errors
E_xkp1_est  = repmat(Xk, 1, n_obs) - cell2mat(Xkp1_est);
E_xkp1_est  = reshape(E_xkp1_est, [], n, n_obs);
mses_xkp1_est = squeeze(mean(E_xkp1_est.^2, [1 2]))';
E_xk_est  = repmat(Xk, 1, n_obs) - cell2mat(Xk_est);
E_xk_est  = reshape(E_xk_est, [], n, n_obs);
mses_xk_est = squeeze(mean(E_xk_est.^2, [1 2]))';

% Calculate output estimation errors and mean-squared errors
E_yk_est  = repmat(Yk, 1, n_obs) - cell2mat(Yk_est);
mses_yk_est = mean(E_yk_est.^2);
E_ykp1_est  = repmat(Yk, 1, n_obs) - cell2mat(Ykp1_est);
mses_yk_est = mean(E_ykp1_est.^2);

% TODO: Need to check these are correct
% Observers: KFSS_old, KFPSS, KFFSS, KF_old, KFP, KFF
assert(isequaln(round(mses_xkp1_est, 6), ...
    [3.590839 3.590839 4.531142 4.531142 4.531142] ...
))
assert(isequaln(round(mses_xk_est, 6), ...
    [nan 2.118641 nan nan 3.106096] ...
))
%TODO: Why are dynamic KFs worse?

% MSE from benchmark simulation:
% >> mse = mean(( ...
%        (bench_sim_results{:, {'xNprocess_1', 'xNprocess_2'}} ...
%         - bench_sim_results{:, {'xNkalman_1', 'xNkalman_2'}} ...
%    )).^2, [1 2])
% 
% mse =
% 
%     3.5908
% 
assert(all(round(mses_xkp1_est(1:2), 6) == 3.590839))


%% SISO system example from GEL-7029 course, homework 12.
% See file /gel-7029/homework/hw12/hw12_p3_kalman.m

sys_test_siso2

% Check if benchmark simulation data file exists
if ~isfile('results/hw12_p3_kalman_sim_benchmark.csv')
    error("Run 'kalman_benchmark_hw12_p3.m' to generate benchmark data.")
end

% Observer parameters
W = 0.5; % estimate of Wp used in the filter design
V = 0.8; % estimate of Vp used in the filter design
P0 = eye(n)*1000; % Initialize covariance matrix
Q = diag(repmat(W,n,1));
R = diag(repmat(V,ny,1));
x0_est = [0.1; 0.5];

% Define dynamic Kalman filter
label = 'KF';
KF_old = KalmanFilter(A,B,C,Ts,P0,Q,R,label,x0_est);
assert(strcmp(KF_old.type, "KF"))
assert(isequal(KF_old.A, A))
assert(isequal(KF_old.B, B))
assert(isequal(KF_old.C, C))
assert(isequal(KF_old.Ts, Ts))
assert(isequal(KF_old.P0, P0))
assert(isequal(KF_old.Q, Q))
assert(isequal(KF_old.R, R))
assert(isequal(KF_old.Pkp1, P0))
assert(isequal(KF_old.label, label))
assert(isequal(KF_old.xkp1_est, x0_est))
assert(KF_old.ykp1_est == C*x0_est)
assert(KF_old.n == n)
assert(KF_old.nu == nu)
assert(KF_old.ny == ny)

% New version
label = 'KFF';
model.A = A;
model.B = B;
model.C = C;
model.Q = Q;
model.R = R;
model.Ts = Ts;
KFF = KalmanFilterF(model,P0,label,x0_est);
assert(strcmp(KFF.type, "KFF"))
assert(isequal(KFF.model, model))
assert(isequal(KFF.label, label))
assert(isequal(KFF.xkp1_est, x0_est))
assert(KFF.ykp1_est == C*x0_est)
assert(isequal(KFF.Pkp1, P0))
assert(all(isnan(KFF.xk_est)))
assert(all(isequaln(KFF.Pk, nan(2))))
assert(all(isnan(KFF.yk_est)))
assert(KFF.n == n)
assert(KFF.nu == nu)
assert(KFF.ny == ny)

% Re-define with no initial state specified (should be set to zero)
KF_old = KalmanFilter(A,B,C,Ts,P0,Q,R);
assert(all(isnan(KF_old.K)))
assert(isequal(KF_old.Pkp1, P0))
assert(isequal(KF_old.xkp1_est, zeros(n, 1)))
assert(KF_old.ykp1_est == 0)

% Re-define with no initial state specified (should be set to zero)
KFF = KalmanFilterF(model,P0);
assert(all(isnan(KFF.Kf)))
assert(isequal(KFF.xkp1_est, zeros(n, 1)))
assert(isequaln(KFF.Pkp1, P0))
assert(KFF.ykp1_est == 0)
assert(all(isnan(KFF.xk_est)))
assert(isequaln(KFF.Pk, nan(2)))
assert(all(isnan(KFF.yk_est)))

% number of points to simulate
nT = 100;

% seed random number generator
rng(0)

% Measurement noise for the whole simulation
V = (sqrt(Rp)*randn(1,nT))';

% Process noise for the whole simulation
W = (sqrt(Qp)*randn(2,nT))';

% Intialize system (at k = 0)
x0_est = zeros(n, 1);
xk = x0_est;

% Input signal - pseudo-random binary sequence
U = nan(nT,1);
warning('off')
U(1:nT,1) = idinput(nT, 'PRBS', [0 0.5]);
warning('on')

% Matrices to collect simulation data
xNprocess = nan(n, nT); % process states
yNprocess = nan(ny, nT); % process outputs
xNkalman1 = nan(n, nT); % estimated states
yNkalman1 = nan(ny, nT); % estimated process outputs
KNkalman1 = nan(n, nT); % observer correction gains
diagPNkalman1 = nan(n, nT); % diag of observer covariance matrix
xNkalman2 = nan(n, nT); % estimated states
yNkalman2 = nan(ny, nT); % estimated process outputs
KNkalman2 = nan(n, nT); % observer correction gains
diagPNkalman2 = nan(n, nT); % diag of observer covariance matrix

t = Ts * (0:nT-1)';

for i = 1:nT

    % Process output in current timestep
    yk = C*xk + V(i,:)';

    % Record process states and output
    xNprocess(:, i) = xk;
    yNprocess(:, i) = yk;

    % Record Kalman filter estimates of current
    % states and process outputs (made in previous
    % timestep)
    xNkalman1(:, i) = KF_old.xkp1_est;
    yNkalman1(:, i) = KF_old.ykp1_est;

    % Check predictions are the same
    assert(all(abs(KF_old.xkp1_est - KFF.xkp1_est) < 1e-12))
    assert(all(abs(KF_old.ykp1_est - KFF.ykp1_est) < 1e-13))

    % Update KFs
    KF_old.update(yk, U(i));
    KFF.update(yk, U(i));

    % Record updated estimates of current states
    % and process outputs (filtering KF only)
    xNkalman2(:, i) = KFF.xk_est;
    yNkalman2(:, i) = KFF.yk_est;

    % Record Kalman filter variables
    KNkalman1(:, i) = KF_old.K;
    diagPNkalman1(:, i) = diag(KF_old.Pkp1);
    KNkalman2(:, i) = KFF.Kf;
    diagPNkalman2(:, i) = diag(KFF.Pkp1);

    % Process states in next timestep
    xk = A*xk + B*U(i) + W(i,:)';

end

% Record final Kalman filter estimates
xNkalman1(:, nT) = KF_old.xkp1_est;
yNkalman1(:, nT) = C * KF_old.xkp1_est;

sim_results = [table(t,U,V,W) ...
    array2table(xNprocess', 'VariableNames', {'x1', 'x2'}) ...
    array2table(xNkalman1', 'VariableNames', {'x1_est_KF', 'x2_est_KF'}) ...
    array2table(yNprocess', 'VariableNames', {'y'}) ...
    array2table(yNkalman1', 'VariableNames', {'y_est_KF'}) ...
    array2table(KNkalman1', 'VariableNames', {'K1', 'K2'}) ...
    array2table(diagPNkalman1', 'VariableNames', {'P1', 'P2'})];

% Verify results by comparing with code from course homework
filename = 'hw12_p3_kalman_sim_benchmark.csv';

warnId = 'MATLAB:table:ModifiedAndSavedVarnames';
warnStruct = warning('off',warnId);
bench_sim_results = readtable(fullfile('results', filename));
warning(warnStruct);

% Compare results to benchmark results
% tail(sim_results)
% tail(bench_sim_results)

% Check states and outputs
assert(isequal( ...
    round(sim_results{1:10, {'x1', 'x2', 'y'}}, 6), ...
    round(bench_sim_results{1:10, {'X_t__1', 'X_t__2', 'y_t_'}}, 6) ...
))

% Check state and output estimates
assert(isequal( ...
    round(sim_results{1:100, {'x1_est_KF', 'x2_est_KF', 'y_est_KF'}}, 6), ...
    round(bench_sim_results{1:100, {'X_e_t__1', 'X_e_t__2', 'y_est_t_'}}, 6) ...
))

% Check correction gains
assert(isequal( ...
    round(sim_results{1:100, {'K1', 'K2'}}, 6), ...
    round(bench_sim_results{1:100, {'K_t__1', 'K_t__2'}}, 6) ...
))

% Check P covariances
assert(isequal( ...
    round(sim_results{1:100, {'P1', 'P2'}}, 6), ...
    round(bench_sim_results{1:100, {'P_t__1', 'P_t__2'}}, 6) ...
))


% plot results

% figure(1); clf
%
% subplot(411);
% plot(t', yNprocess, 'k', t', yNkalman, 'r', 'Linewidth', 2)
% legend('Process output', 'KF estimates')
% ylabel('$y$')
% grid on
%
% subplot(412);
% plot(t', xNprocess(1,:), 'k', t', xNkalman(1,:), 'r', 'Linewidth', 2)
% legend('Process state', 'KF estimates')
% ylabel('$x_1$')
% grid on
%
% subplot(413);
% plot(t', xNprocess(2,:), 'k', t', xNkalman(2,:), 'r', 'Linewidth', 2);
% legend('Process state', 'KF estimates')
% ylabel('$x_2$')
% grid on
%
% subplot(414);
% stairs(t', U', 'Linewidth', 2);
% xlabel('Time [s]');
% ylabel('$u_1$')
% grid on


%% Test on 2x2 system

% Check if benchmark simulation data file exists
filename = 'KF_sim_benchmark_2x2.csv';
if ~isfile(fullfile(results_dir, filename))
    error("Test data file '%s' not found.", filename)
end

bench_sim_results = readtable(fullfile(results_dir, filename));

% Noise variances
sigma_p = 0;  % input disturbances
sigma_M = 0.1;  % measurement noise

% Sample time
Ts = 1;

% Discrete time state space model
A = [ 0.8890       0     1 -0.2;
           0  0.8890  -0.2    1;
           0       0     1    0;
           0       0     0    1];
B = [    1 -0.2  0  0;  % TODO: increase the coupling, -0.5?
      -0.2    1  0  0;
         0    0  1  0;
         0    0  0  1];
C = [ 0.1110 0         0  0;
             0  0.1110 0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Designate known inputs and measured outputs
u_meas = [1 2];
y_meas = [1 2];
np = nu - numel(u_meas);
nu = nu - np;
Bu = B(:, u_meas);
Du = D(:, u_meas);

% Default initial condition
x0_est = zeros(n, 1);

% Observer parameters
Q = 0.01 .* eye(n);
R = 0.1 .* eye(ny);
P0 = 1e-4 .* eye(n);

% Prepare a struct of model parameters
model = struct();
model.A = A;
model.B = Bu;
model.C = C;
model.Ts = Ts;
model.Q = Q;
model.R = R;

% Define Kalman filters
KFPSS = KalmanFilterPSS(model,"KFPSS",x0_est);
KFFSS = KalmanFilterFSS(model,"KFFSS",x0_est);
KFP = KalmanFilterP(model,P0,"KFP",x0_est);
KFF = KalmanFilterF(model,P0,"KFF",x0_est);

KF_old2 = KalmanFilter(A,Bu,C,Ts,P0,Q,R,"KF_old2",x0_est);

observers = {KFPSS, KFFSS, KFP, KFF};
obs_labels = cellfun(@(obs) obs.label, observers);
n_obs = numel(observers);

% Number of timesteps to simulate
nT = 100;
t = Ts*(0:nT-1)';

% Set RNG seed
rng(0)

% Random inputs
P = sigma_p^2 .* [zeros(5, np); randn([45 np]); zeros(nT-50, np)];
U = [zeros(5, nu); idinput([45 nu]); [-0.1 0.25].*ones(nT-50, nu)];
V = sigma_M^2 .* randn([nT+1 ny]);
V = V(1:nT, :);

% Simulate system
U_sim = [U P];
[Yk, t, Xk] = lsim(Gpss, U_sim, t, x0_est);

% Add measurement noise
Yk_m = Yk + V;

% Plot results

% figure(2); clf
%
% ax1 = subplot(311);
% plot(t, Y_m(:, 1), 'o', t, Y_m(:, 2), 'o'); hold on
% % Modify plot colors
% cols = get(gca,'ColorOrder');
% cols(ny+1:2*ny, :) = cols(1:ny, :);
% set(gca,'ColorOrder', cols);
% plot(t, Y(:, 1), t, Y(:, 2), 'Linewidth', 2)
% legend({'$y_{m,1}(k)$', '$y_{m,2}(k)$', '$y_1(k)$', '$y_2(k)$'}, 'Interpreter', 'latex')
% ylabel('$y_i(k)$', 'Interpreter', 'latex')
% title('Outputs')
% grid on
%
% ax2 = subplot(312);
% for i = 1:nu
%     stairs(t, U(:, i), 'Linewidth', 2); hold on
% end
% ylim([-1.2 1.2])
% legend({'$u_1(k)$', '$u_2(k)$'}, 'Interpreter', 'latex')
% ylabel('$u_i(k)$', 'Interpreter', 'latex')
% title('Known inputs')
% grid on
%
% ax3 = subplot(313);
% for i = 1:np
%     stairs(t, P(:, i), 'Linewidth', 2); hold on
% end
% legend({'$p_1(k)$', '$p_2(k)$'}, 'Interpreter', 'latex')
% ylabel('$p_i(k)$', 'Interpreter', 'latex')
% title('Unknown inputs')
% grid on

% Matrices to collect simulation data
Xk_est = cell(1, n_obs);
Yk_est = cell(1, n_obs);
Xkp1_est = cell(1, n_obs);
Ykp1_est = cell(1, n_obs);
trPk_obs = cell(1, n_obs);
trPkp1_obs = cell(1, n_obs);
for f = 1:n_obs
    Xk_est{f} = nan(nT, n);
    Yk_est{f} = nan(nT, ny);
    Xkp1_est{f} = nan(nT, n);
    Ykp1_est{f} = nan(nT, ny);
end

for i = 1:nT

    yk = Yk_m(i, :)';
    if i > 1
        uk = U(i-1, :)';
    else
        uk = zeros(nu, 1);
    end

    for f = 1:n_obs

        obs = observers{f};

        % Record Kalman filter predictions
        if isprop(obs, "xkp1_est") || isfield(obs, "xkp1_est")
            Xkp1_est{f}(i, :) = obs.xkp1_est';
        end
        if isprop(obs, "ykp1_est") || isfield(obs, "ykp1_est")
            Ykp1_est{f}(i, :) = obs.ykp1_est';
            assert(isequaln(C * obs.xkp1_est, obs.ykp1_est))
        end
        switch obs.type
            case {"KFP", "KFPSS"}
                trPkp1_obs{f}(i, :) = trace(obs.Pkp1);
        end

        % Update Kalman filter
        if any(strcmp(obs.type, ["KFSS", "KF"]))
            obs = update_KF(obs, uk, yk);
        else
            obs.update(yk, uk);
        end

        % Record Kalman filter updated estimates
        if isprop(obs, "xk_est") || isfield(obs, "xk_est")
            Xk_est{f}(i, :) = obs.xk_est';
        end
        if isprop(obs, "yk_est") || isfield(obs, "yk_est")
            Yk_est{f}(i, :) = obs.yk_est';
            assert(isequaln(C * obs.xk_est, obs.yk_est))
        end
        switch obs.type
            case {"KFF", "KFFSS"}
                trPk_obs{f}(i, :) = trace(obs.Pkp1);
        end

        % TODO: Remove this once not using structs any more
        observers{f} = obs;

    end

end

% Calculate state estimation errors and mean-squared errors
E_xkp1_est  = repmat(Xk, 1, n_obs) - cell2mat(Xkp1_est);
E_xkp1_est  = reshape(E_xkp1_est, [], n, n_obs);
mses_xkp1_est = squeeze(mean(E_xkp1_est.^2, [1 2]))';
E_xk_est  = repmat(Xk, 1, n_obs) - cell2mat(Xk_est);
E_xk_est  = reshape(E_xk_est, [], n, n_obs);
mses_xk_est = squeeze(mean(E_xk_est.^2, [1 2]))';

% Calculate output estimation errors and mean-squared errors
E_ykp1_est  = repmat(Yk, 1, n_obs) - cell2mat(Ykp1_est);
E_ykp1_est  = reshape(E_ykp1_est, [], ny, n_obs);
mses_ykp1_est = mean(E_ykp1_est.^2);
E_yk_est  = repmat(Yk, 1, n_obs) - cell2mat(Yk_est);
E_yk_est  = reshape(E_yk_est, [], ny, n_obs);
mses_yk_est = mean(E_yk_est.^2);


% Save benchmark results - from observers 3, 4
% bench_sim_results = [table(t,U,P,Y,V,Y_m) ...
%     array2table(X, 'VariableNames', {'X_1', 'X_2', 'X_3', 'X_4'}) ...
%     array2table(Xkp1_est(:, 2*n+1:4*n), 'VariableNames', {'X_est_1', 'X_est_2', 'X_est_3', 'X_est_4', 'X_est_5', 'X_est_6', 'X_est_7', 'X_est_8'}) ...
%     array2table(Ykp1_est(:, 2*ny+1:4*ny), 'VariableNames', {'Y_est_1', 'Y_est_2', 'Y_est_3', 'Y_est_4'}) ...
%     array2table(trP_obs(:, 3:4), 'VariableNames', {'trP_obs_1', 'trP_obs_2'})];
% writetable(bench_sim_results, fullfile(results_dir, filename))

% Combine results into table - only first two observers
% sim_results = [table(t,U,P,Yk,V,Yk_m) ...
%     array2table(Xk, 'VariableNames', {'X_1', 'X_2', 'X_3', 'X_4'}) ...
%     array2table(Xkp1_est(:, 1:2*n), 'VariableNames', {'X_est_1', 'X_est_2', 'X_est_3', 'X_est_4', 'X_est_5', 'X_est_6', 'X_est_7', 'X_est_8'}) ...
%     array2table(Ykp1_est(:, 1:2*ny), 'VariableNames', {'Y_est_1', 'Y_est_2', 'Y_est_3', 'Y_est_4'}) ...
%     array2table(trP_obs(:, 1:2), 'VariableNames', {'trP_obs_1', 'trP_obs_2'})];
%sim_results

% Display results
sim_results = [
    table(t,U,P,Yk,V,Yk_m) ...
    array2table(Xk, 'VariableNames', {'x1', 'x2', 'x3', 'x4'}) ...
];
for f = 1:n_obs
    obs = observers{f};
    labels = compose("xkp1_est_%d_", 1:n) + obs.label;
    sim_results(:, labels) = array2table(Xkp1_est{f});
end
for f = 1:n_obs
    obs = observers{f};
    labels = compose("xk_est_%d_", 1:n) + obs.label;
    sim_results(:, labels) = array2table(Xk_est{f});
end
for f = 1:n_obs
    obs = observers{f};
    labels = compose("yk_est_%d_", 1:ny) + obs.label;
    sim_results(:, labels) = array2table(Yk_est{f});
end
for f = 1:n_obs
    obs = observers{f};
    labels = compose("ykp1_est_%d_", 1:ny) + obs.label;
    sim_results(:, labels) = array2table(Ykp1_est{f});
end
%head(sim_results)

% Plot observer estimates
% j = 4;
% figure(3); clf
% plot(t, Ykp1_est(:, (j-1)*ny+1), 'o', t, Ykp1_est(:, (j-1)*ny+2), 'o'); hold on
% % Modify plot colors
% cols = get(gca,'ColorOrder');
% cols(ny+1:2*ny, :) = cols(1:ny, :);
% set(gca,'ColorOrder', cols);
% plot(t, Y(:, 1), t, Y(:, 2), 'Linewidth', 2)
% xlim(t([1 end]))
% legend({'$\hat{y}_1(k)$', '$\hat{y}_2(k)$', '$y_1(k)$', '$y_2(k)$'}, 'Interpreter', 'latex')
% ylabel('$y_i(k)$', 'Interpreter', 'latex')
% title('Outputs')
% grid on

% Verify results by comparing with saved benchmark results
%head(bench_sim_results)
% First check simulation inputs and outputs
selected_labels = {'t', 'U', 'P', 'Yk', 'V', 'Yk_m', 'x1', 'x2', 'x3', 'x4'};
assert(isequaln( ...
    round(sim_results{:, selected_labels}, 7), ...
    round(bench_sim_results{1:100, 1:15}, 7) ...
))

% Check observer estimates - steady state Kalman filters
selected_labels = ["xkp1_est_1" "xkp1_est_2" "xkp1_est_3" "xkp1_est_4"];
assert(isequaln( ...
    round(sim_results{:, selected_labels + "_" + "KFPSS"}, 7), ...
    round(bench_sim_results{1:100, {'X_est_1', 'X_est_2', 'X_est_3', 'X_est_4'}}, 7) ...
))
assert(isequaln( ...
    round(sim_results{:, selected_labels + "_" + "KFFSS"}, 7), ...
    round(bench_sim_results{1:100, {'X_est_1', 'X_est_2', 'X_est_3', 'X_est_4'}}, 7) ...
))

% Check observer estimates - dynamic Kalman filters
selected_labels = ["xkp1_est_1" "xkp1_est_2" "xkp1_est_3" "xkp1_est_4"];
assert(isequaln( ...
    round(sim_results{:, selected_labels + "_" + "KFP"}, 7), ...
    round(bench_sim_results{1:100, {'X_est_5', 'X_est_6', 'X_est_7', 'X_est_8'}}, 7) ...
))
assert(isequaln( ...
    round(sim_results{:, selected_labels + "_" + "KFF"}, 7), ...
    round(bench_sim_results{1:100, {'X_est_5', 'X_est_6', 'X_est_7', 'X_est_8'}}, 7) ...
))


%% Test copy methods

sys_test_siso

% Covariance matrices
Q = diag([0.1; 0.1]);
R = 0.5;
N = zeros(n,ny);
x0_est = [0.1; 0.5];
P0 = diag([1e-4 1e-4]);

% Prepare a struct of model parameters
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.Ts = Ts;
model.Q = Q;
model.R = R;

% Define steady-state Kalman
KFFSS = KalmanFilterPSS(model,"KFSS",x0_est);

% Define dynamic Kalman filter
KF_old = KalmanFilter(A,B,C,Ts,P0,Q,R,"KF");

% Define dynamic Kalman filter - filtering form
model.A = A;
model.B = B;
model.C = C;
model.Q = Q;
model.R = R;
model.Ts = Ts;
KFF = KalmanFilterF(model,P0,"KFF");

% Test handle copy
KFSS_hcopy = KFFSS;
assert(isequaln(KFSS_hcopy, KFFSS))  % same values
assert(KFSS_hcopy == KFFSS)  % must be same object

KFFSS.x0 = [0.2; 0.5];
assert(isequal(KFSS_hcopy.x0, [0.2; 0.5]))

KF_hcopy = KF_old;
assert(isequaln(KF_hcopy, KF_old))  % same values
assert(KF_hcopy == KF_old)  % must be same object

KF_old.label = "New name";
assert(isequal(KF_hcopy.label, "New name"))

KFF_hcopy = KFF;
assert(isequaln(KFF_hcopy, KFF))  % same values
assert(KFF_hcopy == KFF)  % must be same object

KFF.model.A(1, 1) = KFF.model.A(1, 1) + 0.1;
assert(isequal(KFF_hcopy.model.A(1, 1), KFF.model.A(1, 1)))

% Redefine steady-state Kalman
KFFSS = KalmanFilterPSS(model,"KFSS",x0_est);

% Redefine dynamic Kalman filter
KF_old = KalmanFilter(A,B,C,Ts,P0,Q,R,"KF");

% Redefine dynamic Kalman filter
KFF = KalmanFilterF(model,P0,"KFF");

% Test true copy
KFSS_copy = KFFSS.copy();
assert(isequaln(KFSS_copy, KFFSS))  % same values
assert(KFSS_copy ~= KFFSS)  % must not be same object

KFFSS.x0 = [0.2; 0.5];
assert(~isequal(KFSS_copy.x0, [0.2; 0.5]))

KF_copy = KF_old.copy();
assert(isequaln(KF_copy, KF_old))  % same values
assert(KF_copy ~= KF_old)  % must not be same object

KF_old.label = "New name";
assert(~isequal(KF_copy.label, "New name"))

KFF_copy = KFF.copy();
assert(isequaln(KFF_copy, KFF))  % same values
assert(KFF_copy ~= KFF)  % must not be same object

KFF.model.A(1, 1) = KFF.model.A(1, 1) + 0.1;
assert(~isequal(KFF_copy.model.A, KFF.model.A(1, 1)))
