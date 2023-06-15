% Test MKFObserverSF_RODD95 class

clear all

% This is only needed for plotting
addpath("../plot-utils/")
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir)
end

seed = 0;
rng(seed)


%% Test sequence generation 1

% Load SISO system and disturbance model from file
sys_rodin_step

% Define sequence fusion observer
% This example is used in methods section of thesis report
P0 = eye(n);
Q0 = diag([0.01 0]);
R = sigma_M^2;
f = 9;  % fusion horizon
m = 1;  % maximum number of shocks
d = 3;  % spacing parameter
io = struct('u_known', u_known, 'y_meas', true(ny,1));
MKF_SF95 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,"MKF_SF_RODD95");

seq_test = { ...
    [1 1 1 1 1 1 1 1 1];
    [2 1 1 1 1 1 1 1 1];
    [1 1 1 2 1 1 1 1 1];
    [1 1 1 1 1 1 2 1 1]
};
assert(isequal(MKF_SF95.seq, seq_test))

% Check hypothesis probabilities
alpha = 1 - [1-epsilon].^d;
assert(MKF_SF95.alpha == alpha)
assert(isequal(MKF_SF95.p_rk, [(1-alpha) alpha]'))
p_seq = [(1-alpha).^3 alpha*(1-alpha).^2 alpha*(1-alpha).^2 alpha*(1-alpha).^2]';
assert(isequal(MKF_SF95.p_seq, p_seq))
assert(MKF_SF95.beta == sum(p_seq));


%% Test sequence generation 2

% Load SISO system and disturbance model from file
sys_rodin_step

% Define sequence fusion observer
P0 = eye(n);
Q0 = diag([0.01 0]);
R = sigma_M^2;
f = 10;  % fusion horizon
m = 2;  % maximum number of shocks
d = 5;  % spacing parameter
io = struct('u_known', u_known, 'y_meas', true(ny,1));
MKF_SF95 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,"MKF_SF_RODD95");

seq_test = { ...
    [1 1 1 1 1 1 1 1 1 1];
    [2 1 1 1 1 1 1 1 1 1];
    [1 1 1 1 1 2 1 1 1 1];
    [2 1 1 1 1 2 1 1 1 1]
};
assert(isequal(MKF_SF95.seq, seq_test))

% Check hypothesis probabilities
alpha = 1 - [1-epsilon].^d;
assert(MKF_SF95.alpha == alpha)
assert(isequal(MKF_SF95.p_rk, [(1-alpha) alpha]'))
p_seq = [(1-alpha).^2 alpha*(1-alpha) alpha*(1-alpha) alpha.^2]';
assert(max(abs(MKF_SF95.p_seq - p_seq)) < 1e-15)
assert(MKF_SF95.beta == sum(p_seq));


%% Test initialization with SISO system

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Check observer attributes
assert(strcmp(MKF_SF95.type, "MKF_SF_RODD95"))
assert(MKF_SF95.epsilon == epsilon)
assert(isequal(MKF_SF95.sigma_wp, sigma_wp))
assert(MKF_SF95.f == 15)
assert(MKF_SF95.m == 1)
assert(MKF_SF95.d == 3)
assert(MKF_SF95.nf == f)
assert(MKF_SF95.nm == 6)
assert(MKF_SF95.nh == 8)
assert(MKF_SF95.nh_max == 8)
assert(size(MKF_SF95.filters.Xkp1_est, 3) == MKF_SF95.nh)
assert(isequaln(MKF_SF95.i, 15))
assert(MKF_SF95.n == 2)
assert(MKF_SF95.nu == 1)
assert(MKF_SF95.ny == 1)
assert(MKF_SF95.nj == 2)
assert(isequal(MKF_SF95.models{1}.A, A) && isequal(MKF_SF95.models{2}.A, A))
assert(isequal(MKF_SF95.models{1}.B, Bu) && isequal(MKF_SF95.models{2}.B, Bu))
assert(isequal(MKF_SF95.models{1}.C, C) && isequal(MKF_SF95.models{2}.C, C))
assert(isequal(MKF_SF95.models{1}.D, Du) && isequal(MKF_SF95.models{2}.D, Du))
assert(MKF_SF95.Ts == Ts)
assert(isequaln(MKF_SF95.io.u_known, u_known))
assert(isequal(MKF_SF95.Q0, Q0))
assert(isequal(MKF_SF95.models{1}.Q, [0.01 0; 0 sigma_wp{1}(1)^2]))
assert(isequal(MKF_SF95.models{2}.Q, [0.01 0; 0 sigma_wp{1}(2)^2]))
assert(isequal(MKF_SF95.models{1}.R, R) && isequal(MKF_SF95.models{2}.R, R))
assert(isequal(size(MKF_SF95.seq), [MKF_SF95.nm 1]))
assert(isequal(MKF_SF95.seq, {
    [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
    [2 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
    [1 1 1 2 1 1 1 1 1 1 1 1 1 1 1]
    [1 1 1 1 1 1 2 1 1 1 1 1 1 1 1]
    [1 1 1 1 1 1 1 1 1 2 1 1 1 1 1]
    [1 1 1 1 1 1 1 1 1 1 1 1 2 1 1] ...
}))
assert(isequal(size(cell2mat(MKF_SF95.seq)), [MKF_SF95.nm MKF_SF95.nf]))
assert(MKF_SF95.beta == sum(MKF_SF95.p_seq))
assert(MKF_SF95.nf == size(MKF_SF95.seq{1}, 2))
assert(all(abs(MKF_SF95.filters.Xkp1_est(:, :) - zeros(n, 1)) == 0, [1 2]))
assert(isequal(MKF_SF95.P0, P0))
assert(isequal(MKF_SF95.epsilon, 0.01))
assert(isequal(round(MKF_SF95.alpha, 6), 0.029701))
assert(isequal(round(MKF_SF95.p_rk, 6), [0.970299 0.029701]'))

% Check optional definition with an initial state estimate
label = 'MKF_SF95_testx0';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 9;  % fusion horizon
m = 2;  % maximum number of shocks
d = 3;  % spacing parameter
x0 = [0.1; 0.5];
MKF_testx0 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0);
assert(all(abs(MKF_testx0.filters.Xkp1_est(:, :) - x0) == 0, [1 2]))
assert(isequal(MKF_testx0.seq, {
    [1 1 1 1 1 1 1 1 1]
    [2 1 1 1 1 1 1 1 1]
    [1 1 1 2 1 1 1 1 1]
    [1 1 1 1 1 1 2 1 1]
    [2 1 1 2 1 1 1 1 1]
    [2 1 1 1 1 1 2 1 1]
    [1 1 1 2 1 1 2 1 1] ...
}))
assert(MKF_testx0.nh == 13)
assert(MKF_testx0.nm == 7)
% Also with initial modes, r(-1)
r0 = 1;
MKF_testx0r0 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0,r0);
assert(all(abs(MKF_testx0r0.filters.Xkp1_est(:, :) - x0) == 0, [1 2]))
assert(isequal(MKF_testx0r0.r0, ones(13, 1)))
r0 = [1 2 1 2 1 2 1 2 1 2 1 2 1]';
MKF_testx0r0 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0,r0);
assert(all(abs(MKF_testx0r0.filters.Xkp1_est(:, :) - x0) == 0, [1 2]))
assert(isequal(MKF_testx0r0.r0, r0))


%% Test observers on SISO system with 1 shock

% This test is a copy of that used for the test simulations in
% thesis report but without plot generation. See this script for
% the latest version:
%  - disturb-models/robertson-1995/rod_obs_test_sim.m

% Load system and disturbance model from file
sys_rodin_step

% Simulation settings
nT = 100;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = 9.5;
du0 = 1;
% When you make the shock larger the MKF observers
% do better
%du0 = 2;

% Measured input
%U = (idinput(size(t)) + 1)/2;
U = zeros(size(t));
U(t >= 1) = -1;

% Disturbance input
alpha = zeros(nT+1, 1);
alpha(t == t_shock) = 1;  % this is used by the SKF observer
Wp = du0 .* alpha;

% Observer model without disturbance noise input
Bu = B(:, u_known);
Du = D(:, u_known);
nu = sum(u_known);
nw = sum(~u_known);

% Set noise variances for observer design
sigma_M = 0.1;
sigma_W = [0; 0];

% Load observers from file
% Various observers designed for the system defined in:
%  - sys_rodin_step.m

addpath("../process-observers")

% Check observability of system
Qobs = obsv(A, C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Observer model without disturbance noise input
Bu = B(:, u_known);
Bw = B(:, ~u_known);
Du = D(:, u_known);

% Specify covariance for state variable 1
% This is used by all observers
q1 = 0.01;

% Multiple model observer with sequence fusion #1
label = 'MKF_SF1';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 15;  % fusion horizon
m = 1;  % maximum number of shocks
d = 5;  % spacing parameter
io = struct('u_known', u_known, 'y_meas', true(ny,1));
MKF_SF1 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model observer with sequence fusion #2
label = 'MKF_SF2';
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 10;  % fusion horizon
m = 1;  % maximum number of shocks
d = 1;  % spacing parameter
io = struct('u_known', u_known, 'y_meas', true(ny,1));
MKF_SF2 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

U_sim = [U Wp];

% Custom MKF test observer

% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t == t_shock)

% Multiple model filter - two sequences, one empty, one correct
obs_model = struct;
obs_model.A = A;
obs_model.B = Bu;
obs_model.C = C;
obs_model.D = Du;
obs_model.Ts = Ts;
obs_models = {obs_model, obs_model};
P0 = 1000*eye(n);
Q0 = diag([0.01 1]);
%P0_init = repmat({P0}, 1, 2);
obs_models{1}.Q = diag([Q0(1,1) sigma_wp{1}(1)^2]);
obs_models{2}.Q = diag([Q0(1,1) sigma_wp{1}(2)^2]);
obs_models{1}.R = sigma_M.^2;
obs_models{2}.R = sigma_M.^2;
seq = {ones(1, nT+1); ones(1, nT+1)};
seq{2}(t == t_shock) = 2;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
MKF3 = MKFObserverS(obs_models,P0,seq,T,'MKF3');

% Multiple model filter - one sequence with correct shock
seq = {ones(1, nT+1)};
seq{1}(t == t_shock) = 2;
p_gamma = [1-epsilon epsilon]';
T = repmat(p_gamma', 2, 1);
MKF4 = MKFObserverS(obs_models,P0,seq,T,'MKF4');

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = SKFObserverS(obs_models,P0,seq{1},"SKF");

% Choose observers to test - plots results of last one
observers = {SKF, MKF3, MKF4, MKF_SF2, MKF_SF1};
%observers = {MKF_SF1, MKF_SF2};

% Note: KF1 is too slow to pass static error test here

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
sigma_MP = 0;  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(size(Y));

% Simulate observers

% Measured inputs (not including disturbances)
U_m = U;

n_obs = numel(observers);
MSE = containers.Map();
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_results.X_est(1:20,:) - X(1:20, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_results.Y_est(1:20,:) - Y(1:20, :)) < 1e-10))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
    if all(sigma_MP == 0)
        assert(abs(sim_results.Y_est(end, :) - Y(end, :)) < 1e-3);
        assert(abs(sim_results.X_est(end, 2) - du0) < 1e-3);
    end

    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE(obs.label) = mean((Y_est - Y).^2);

    % Save updated observer
    observers{i} = obs;

    % Save simulation results
    % TODO: Save data for making plots for thesis report
    % sim_results

end

% % Display results of last simulation
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% trP_obs = sim_results.trP_obs;
% table(t,alpha,U,Wp,X,Y,Y_m,X_est,Y_est,E_obs,trP_obs)
% if isfield(sim_results, "trP_obs_f")
%     trP_obs_f = sim_results.trP_obs_f;
% end
% if isfield(sim_results, "K_obs_f")
%     K_obs_f = sim_results.K_obs_f;
% end

% Display gains and trace of covariance matrices of each filter
% table(t, cell2mat(K_obs_f), trP_obs_f, ...
%     'VariableNames',{'t', 'K_1, K_2', 'trace(P)'});

% Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'})

% % Plot of inputs and outputs
% 
% obs_label = obs.label;
% 
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(groot,'defaulttextinterpreter','latex');
% set(groot,'defaultLegendInterpreter','latex');
% 
% figure(1); clf
% colors = get(gca,'colororder');
% ax1 = subplot(4,1,1);
% stairs(t,Y_m); hold on
% stairs(t,Y_est,'Linewidth',2);
% ax1.ColorOrder = colors(1:size(Y_m,2),:);
% max_min = [min(min([Y_m Y_est])) max(max([Y_m Y_est]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$y_m(k)$ and $\hat{y}(k)$')
% title(sprintf('%s - Process output measurements and estimates', ...
%     escape_latex_chars(obs_label)))
% legend('$y_m(k)$','$\hat{y}(k)$')
% grid on
% 
% ax2 = subplot(4,1,2);
% stairs(t,X); hold on
% stairs(t,X_est,'Linewidth',2);
% ax2.ColorOrder = colors(size(Y,2)+1:size(Y,2)+size(X,2),:);
% max_min = [min(min([X X_est])) max(max([X X_est]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$x_i(k)$ and $\hat{x}_i(k)$')
% labels = repmat({''}, 1, n*2);
% for i=1:n
%     labels{i} = sprintf("$x_{%d}(k)$", i);
% end
% for i=1:n
%     labels{i+n} = sprintf("$%s{x}_{%d}(k)$", '\hat', i);
% end
% legend(labels)
% title('States and state estimates')
% grid on
% 
% ax3 = subplot(4,1,3);
% stairs(t,U,'Linewidth',2); hold on;
% stairs(t,Wp,'Linewidth',2)
% max_min = [min(min([U Wp])) max(max([U Wp]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$u(k)$ and $w_p(k)$')
% legend('$u(k)$', '$w_p(k)$')
% title('Process inputs')
% grid on
% 
% ax4 = subplot(4,1,4);
% stairs(t,alpha,'Color',colors(end,:),'Linewidth',2)
% max_min = [min(min(alpha)) max(max(alpha))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% xlabel('t')
% ylabel('$\gamma(k)$')
% title('Random shock sequence')
% grid on
% 
% linkaxes([ax1, ax2 ax3 ax4], 'x')
% 
% set(gcf,'Position',[100 200 560 600]);


% % Plot of conditional filter probabilities
% switch obs.type
%     case {"MKF_S", "MKF_BM", "MKF_GPB1", "MKF_GPB2", 
%           "MKF_SF", "MKF_SF_DI", "MKF_SF_RODD", "MKF_SF_RODD95"}
%         p_seq_g_Yk = sim_results.MKF_p_seq_g_Yk;
% 
%         figure(11); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Pr(\Gamma(k) \mid Y(k))$'};
%         make_waterfall_plot(t, p_seq_g_Yk, [0 1], ...
%             ax_labels, [0 72]);
%         filename = sprintf('rod_MKFObserver_test_pyk_wfplot');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Conditional probabilities of y(k)')
% end
% 
% % Plot of trace of filter covariance matrices
% switch obs.type
%     case {"MKF_S", "MKF_BM", "MKF_GPB1", "MKF_GPB2", 
%           "MKF_SF", "MKF_SF_DI", "MKF_SF_RODD", "MKF_SF_RODD95"}
%         trP_obs_f = sim_results.trP_obs_f;
% 
%         figure(12); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Tr(P(k))$'};
%         make_waterfall_plot(t, trP_obs_f, [0 2], ax_labels, [0 82]);
%         filename = sprintf('rod_MKFObserver_test_trP_wfplot');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Trace of covariance matrices')
% 
% end
% 
% % Plot of filter correction gains (k1)
% 
% switch obs.type
%     case {'MKF_SF'}
%         K_obs_f = cell2mat(sim_results.K_obs_j);
%         % Select first gain value onlu
%         K1_obs = K_obs_f(:,1:2:end);
%         K2_obs = K_obs_f(:,2:2:end);
% 
%         figure(13); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$K_{1,1}(k)$'};
%         make_waterfall_plot(t, K1_obs, [0 5], ax_labels, [0 82]);
%         title('Filter correction gain ($K_{1,1}$)', 'Interpreter', 'latex')
%         filename = sprintf('rod_MKFObserver_test_K_wfplot');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
% 
% end

% Test MSE values
% Results on 2022-11-29 after fixing MKF_SF_DI and MKF_SF_RODD95
MSE_test_values = containers.Map(...
    {'KF2',   'KF3',   'MKF_SF1',  'MKF_SF2',  'MKF3',  'MKF4',  'SKF'}, ...
    [0.000006 0.000992 0.000729 0.000489 0.000499 0.000012 0.000012] ...
);
% for label = MSE.keys
%     fprintf("%s: %f (%f)\n", label{1}, MSE(label{1}), MSE_test_values(label{1}))
% end
for label = MSE.keys
    assert(isequal(round(MSE(label{1}), 6), MSE_test_values(label{1})))
end


%% Test MKF observers on 2x2 system

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

% Designate measured input and output signals
u_known = [true; true; false; false];
y_meas = [true; true];

% Observer model without disturbance noise input
Bu = B(:, u_known);
Du = D(:, u_known);
nu = sum(u_known);
nw = sum(~u_known);
obs_model = struct();
obs_model.A = A;
obs_model.B = Bu;
obs_model.C = C;
obs_model.D = Du;
obs_model.Ts = Ts;

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_known);
nw = sum(~u_known);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.1; 0.1];
sigma_wp = {[0.01 1], [0.01 1]};

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
P0 = 1000*eye(n);
obs_model1 = obs_model;
obs_model1.Q = diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(1)^2]);
obs_model1.R = diag(sigma_M.^2);
KF1 = KalmanFilterF(obs_model1,P0,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
obs_model2 = obs_model;
obs_model2.Q = diag([0.01 0.01 sigma_wp{1}(2)^2 sigma_wp{2}(2)^2]);
obs_model2.R = diag(sigma_M.^2);
KF2 = KalmanFilterF(obs_model2,P0,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
obs_model3 = obs_model;
obs_model3.Q = diag([0.01 0.01 0.1^2 0.1^2]);
obs_model3.R = diag(sigma_M.^2);
KF3 = KalmanFilterF(obs_model3,P0,'KF3');

% System model parameter struct
sys_model = struct();
sys_model.A = A;
sys_model.B = B;
sys_model.C = C;
sys_model.D = D;
sys_model.Ts = Ts;

% Input output signals
io = struct;
io.u_known = u_known;
io.y_meas = y_meas;

% Multiple model filter 1
label = 'MKF_SF1';
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 6;  % fusion horizon
m = 1;  % maximum number of shocks
d = 2;  % spacing parameter
MKF_SF1 = MKFObserverSF_RODD95(sys_model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model filter 2
label = 'MKF_SF2';
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 10;  % fusion horizon
m = 2;  % maximum number of shocks
d = 2;  % spacing parameter
MKF_SF2 = MKFObserverSF_RODD95(sys_model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Check observer initialization
assert(isequal(MKF_SF1.epsilon, epsilon))
assert(isequal(MKF_SF1.sigma_wp, sigma_wp))
assert(MKF_SF1.nm == 7)
assert(MKF_SF1.nh == 13)
assert(isequaln(MKF_SF1.i, 6))
assert(MKF_SF1.n == 4)
assert(MKF_SF1.nu == 2)
assert(MKF_SF1.ny == 2)
assert(MKF_SF1.nj == 3)
assert(isequal(MKF_SF1.models{1}.A, A) && isequal(MKF_SF1.models{2}.A, A))
assert(isequal(MKF_SF1.models{1}.B, Bu) && isequal(MKF_SF1.models{2}.B, Bu))
assert(isequal(MKF_SF1.models{1}.C, C) && isequal(MKF_SF1.models{2}.C, C))
assert(isequal(MKF_SF1.models{1}.D, Du) && isequal(MKF_SF1.models{2}.D, Du))
assert(MKF_SF1.Ts == Ts)
assert(MKF_SF1.models{1}.Ts == Ts)
assert(MKF_SF1.models{2}.Ts == Ts)
assert(isequaln(MKF_SF1.io.u_known, u_known))
assert(isequal(MKF_SF1.Q0, Q0))
assert(isequal(MKF_SF1.models{1}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(1)^2])))
assert(isequal(MKF_SF1.models{2}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(2)^2])))
assert(isequal(MKF_SF1.models{3}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(2)^2 sigma_wp{2}(1)^2])))
assert(isequal(MKF_SF1.R, R))
assert(isequal(MKF_SF1.models{1}.R, R) && ...
    isequal(MKF_SF1.models{2}.R, R) && ...
    isequal(MKF_SF1.models{3}.R, R))
assert(isequal(size(MKF_SF1.seq), [MKF_SF1.nm 1]))
assert(isequal(size(cell2mat(MKF_SF1.seq)), [MKF_SF1.nm MKF_SF1.nf]))
assert(MKF_SF1.beta == sum(MKF_SF1.p_seq))
assert(MKF_SF1.nf == size(MKF_SF1.seq{1}, 2))
assert(isequal(MKF_SF1.P0, P0))
assert(isequal(round(MKF_SF1.alpha, 6), [0.019900 0.019900]'))
p_gamma = [1-MKF_SF1.alpha'; MKF_SF1.alpha'];
Z = [0 0; 0 1; 1 0];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
assert(isequal(round(p_gamma, 6), [0.960977 0.019512 0.019512]'))
assert(isequal(round(MKF_SF1.p_rk, 6), [0.960977 0.019512 0.019512]'))

% Check observer initialization
assert(isequal(MKF_SF2.epsilon, epsilon))
assert(isequal(MKF_SF2.sigma_wp, sigma_wp))
assert(MKF_SF2.nm == 56)
assert(MKF_SF2.nh == 116)
assert(isequaln(MKF_SF2.i, 10))
assert(MKF_SF2.n == 4)
assert(MKF_SF2.nu == 2)
assert(MKF_SF2.ny == 2)
assert(MKF_SF2.nj == 4)
assert(isequal(MKF_SF2.models{1}.A, A) && isequal(MKF_SF2.models{2}.A, A))
assert(isequal(MKF_SF2.models{1}.B, Bu) && isequal(MKF_SF2.models{2}.B, Bu))
assert(isequal(MKF_SF2.models{1}.C, C) && isequal(MKF_SF2.models{2}.C, C))
assert(isequal(MKF_SF2.models{1}.D, Du) && isequal(MKF_SF2.models{2}.D, Du))
assert(MKF_SF2.Ts == Ts)
assert(MKF_SF2.models{1}.Ts == Ts)
assert(MKF_SF2.models{2}.Ts == Ts)
assert(isequaln(MKF_SF2.io.u_known, u_known))
assert(isequal(MKF_SF2.Q0, Q0))
assert(isequal(MKF_SF2.models{1}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(1)^2])))
assert(isequal(MKF_SF2.models{2}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(2)^2])))
assert(isequal(MKF_SF2.models{3}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(2)^2 sigma_wp{2}(1)^2])))
assert(isequal(MKF_SF2.models{4}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(2)^2 sigma_wp{2}(2)^2])))
assert(isequal(MKF_SF2.R, R))
assert(isequal(MKF_SF2.models{1}.R, R) && isequal(MKF_SF2.models{2}.R, R))
assert(isequal(MKF_SF2.models{3}.R, R) && isequal(MKF_SF2.models{4}.R, R))
assert(isequal(size(MKF_SF2.seq), [MKF_SF2.nm 1]))
assert(isequal(size(cell2mat(MKF_SF2.seq)), [MKF_SF2.nm MKF_SF2.nf]))
assert(MKF_SF2.beta == sum(MKF_SF2.p_seq))
assert(MKF_SF2.nf == size(MKF_SF2.seq{1}, 2))
assert(isequal(MKF_SF2.P0, P0))
assert(sum(MKF_SF2.p_rk) == 1)
assert(isequal(round(MKF_SF2.alpha, 6), [0.019900 0.019900]'))
p_gamma = [1-MKF_SF2.alpha'; MKF_SF2.alpha'];
Z = [0 0; 0 1; 1 0; 1 1];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
assert(isequal(round(p_gamma, 6), ...
    [0.960596 0.019504 0.019504 0.000396]'))
assert(isequal(round(MKF_SF2.p_rk, 6), ...
    [0.960596 0.019504 0.019504 0.000396]'))

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5; -0.2; -0.4];
label = 'MKF_testx0';
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
%R = diag([1; 2.3].*sigma_M.^2);
f = 30;  % fusion horizon
m = 2;  % maximum number of shocks
d = 10;  % spacing parameter
MKF_testx0 = MKFObserverSF_RODD95(sys_model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0);
assert(all(abs(MKF_testx0.filters.Xkp1_est(:, :) - x0) == 0, [1 2]))

% Simulation settings
nT = 200;
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
alpha = zeros(nT+1, 2);
alpha(t == t_shock(1), 1) = 1;
alpha(t == t_shock(2), 2) = 1;
Wp = du0' .* alpha;

U_sim = [U Wp];

% Custom MKF test observer
% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)
% Multiple model filter 1
obs_model = struct;
obs_model.A = A;
obs_model.B = Bu;
obs_model.C = C;
obs_model.D = Du;
obs_model.Ts = Ts;
obs_models = {obs_model, obs_model, obs_model};
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 1 1]);
obs_models{1}.Q = diag([Q0(1,1) Q0(2,2) sigma_wp{1}(1)^2 sigma_wp{2}(1)^2]);
obs_models{2}.Q = diag([Q0(1,1) Q0(2,2) sigma_wp{1}(2)^2 sigma_wp{2}(1)^2]);
obs_models{3}.Q = diag([Q0(1,1) Q0(2,2) sigma_wp{1}(1)^2 sigma_wp{2}(2)^2]);
obs_models{1}.R = diag(sigma_M.^2);
obs_models{2}.R = diag(sigma_M.^2);
obs_models{3}.R = diag(sigma_M.^2);
seq = {ones(1, nT+1); ones(1, nT+1); ones(1, nT+1)};
seq{2}(t == t_shock(1)) = 2;
seq{3}(t == t_shock(2)) = 3;
p_gamma = [1-epsilon epsilon]';
Z = [0 0; 0 1; 1 0];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
T = repmat(p_gamma', 3, 1);
MKF3 = MKFObserverS(obs_models,P0,seq,T,'MKF3');

seq = {ones(1, nT+1)};
seq{1}(t == t_shock(1)) = 2;
seq{1}(t == t_shock(2)) = 3;
MKF4 = MKFObserverS(obs_models,P0,seq,T,'MKF4');

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = SKFObserverS(obs_models,P0,seq{1},"SKF");

% Choose observers to test
observers = {KF3, SKF, MKF_SF1, MKF_SF2, MKF3, MKF4};

% Note: KF1 is too slow to pass static error test here

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
Y_m = Y + sigma_MP'.*randn(nT+1, ny);

% Simulate observers

% Measured inputs (not including disturbances)
U_m = U;

n_obs = numel(observers);
MSE = containers.Map();
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_results.X_est(1:5,:) - X(1:5, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_results.Y_est(1:5,:) - Y(1:5, :)) < 1e-10, [1 2]))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
    if all(sigma_MP == 0)
        assert(all(abs(sim_results.Y_est(end, :) - Y(end, :)) < 1e-3, [1 2]));
        assert(all(abs(sim_results.X_est(end, 3:4) - du0) < 1e-3, [1 2]));
    end

    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE(obs.label) = mean((Y_est - Y).^2);
    %fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

    % Save updated observer
    observers{i} = obs;

end


% Display results of last simulation
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% trP_obs = sim_results.trP_obs;
% K_obs_j = sim_results.K_obs_j;
% trP_obs_j = sim_results.trP_obs_j;
% table(t,alpha,U,Wp,X,Y,Y_m,X_est,Y_est,E_obs)

% Display gains and trace of covariance matrix
%K_data = cell2mat(cellfun(@(X) X(:)', K_obs, 'UniformOutput', false));
%table(t, K_data, cell2mat(trP_obs), ...
%    'VariableNames', {'t', 'K{1}, K{2}', 'trace(P{1}), trace(P{2})'})

% Show table of mean-squared errors
% table(MSE.keys', cell2mat(MSE.values'), ...
%     'VariableNames', {'Observer', 'MSE'});

% Results on 2022-11-29 after fixing MKF_SF_DI and MKF_SF_RODD95
MSE_test_values = containers.Map(...
 {'KF3',               'MKF_SF1',           'MKF_SF2', ...
  'MKF3',              'MKF4',              'SKF'}, ...
 {[0.000296 0.000433], [0.000351 0.000374], [0.000261 0.000287], ...
  [0.000971 0.000916], [0.000017 0.000022], [0.000017 0.000022]} ...
);

% for label = MSE.keys
%     fprintf("%s: %f, %f (%f, %f)\n", label{1}, MSE(label{1}), MSE_test_values(label{1}))
% end
for label = MSE.keys
    %fprintf("%s: %f, %f (%f, %f)\n", label{1}, MSE(label{1}), MSE_test_values(label{1}))
    assert(isequal(round(MSE(label{1}), 6), MSE_test_values(label{1})))
end
