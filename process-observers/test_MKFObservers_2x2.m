% Tests the multi-model observers on a 2x2 system with RODD disturbance
% 
%  - KalmanFilterF
%      Single Kalman Filter (used for comparison).
%  - SKFObserver
%      KF with switching system model.
%  - MKFObserver
%      Multiple-model KF with switching system models.
%


clear all

addpath("../plot-utils/")

seed = 0;
rng(seed)


%% Simulation test on 2x2 linear system with RODD disturbances

% Sample time
Ts = 1;

% NOTE: this is a previous version of the system with lower
% coupling (-0.2) and epsilon = [0.01; 0.01].

% Discrete time state space model
A = [
      0.8890       0     1 -0.2
           0  0.8890  -0.2    1
           0       0     1    0
           0       0     0    1
];
B = [
         1 -0.2  0  0;
      -0.2    1  0  0;
         0    0  1  0;
         0    0  0  1
];
C = [ 0.1110 0         0  0;
             0  0.1110 0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
[n, nu, ny] = check_dimensions(A, B, C, D);

% Simulation settings
nT = 50;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = [5 15];
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
%P = cumsum(Wp);
U_sim = [U Wp];
x0 = zeros(n, 1);

% System models
sys_model.A = A;
sys_model.B = B;  % input disturbances unmeasured
sys_model.C = C;
sys_model.Ts = Ts;
sys_models = repmat({sys_model}, 1, 3);

% Choose measurement noise for plant
sigma_MP = [0; 0];  % See below for simulation with noise
V = sigma_MP'.*randn(nT+1, ny);

% Run simulation
[X, Y, Ym] = run_simulation_sys(sys_models,U_sim,V,alpha,nT,x0);

% % Simulate system
% X2 = zeros(nT+1,n);
% Y2 = zeros(nT+1,ny);
% xk = x0;
% for i = 1:nT+1
% 
%     % Inputs
%     uk = U_sim(i,:)';
% 
%     % Compute y(k)
%     yk = C * xk + D * uk;
% 
%     % Store results
%     X2(i, :) = xk';
%     Y2(i, :) = yk';
% 
%     % Compute x(k+1)
%     xk = A * xk + B * uk;
% 
% end
% 
% % Check simulation output is correct
% [Y3, t, X3] = lsim(Gpss, U_sim, t, x0);
% assert(isequal(X, X2))
% assert(isequal(Y, Y2))
% assert(isequal(X, X3))
% assert(isequal(Y, Y3))

% % Plot of inputs and outputs
% figure(4); clf
% P = cumsum(Wp);
% u_labels = [compose("$u_%d(k)$", 1:size(U,2)) compose("$p_%d(k)$", 1:size(P,2))];
% y_labels = [compose("$y_%d(k)$", 1:ny) compose("$y_{m,%d}(k)$", 1:ny)];
% x_label = "Time, $t$";
% make_iodplot(Y,Ym,t,[U P],u_labels,y_labels,x_label)

% Designate measured input and output signals
u_known = [true; true; false; false];
y_meas = [true; true];

% Observer model (without unmeasured disturbance input)
Bu = B(:, u_known);
Du = D(:, u_known);
nu = sum(u_known);
nw = sum(~u_known);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_known);
nw = sum(~u_known);

% Dimensions of observer model
[n, nu, ny] = check_dimensions(A, Bu, C, Du);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.05; 0.05];
%sigma_M = [0; 0];  % set to zero for testing
sigma_wp = [0.01 1; 0.01 1];

% Observer models
model.A = A;
model.B = Bu;  % input disturbances unmeasured
model.C = C;
model.Ts = Ts;
models = repmat({model}, 1, 3);

% Observer parameters (same for all observers)
models{1}.Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
models{2}.Q = diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]);
models{3}.Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2]);
assert(isequal(size(models{1}.Q), [n n]))
assert(isequal(size(models{2}.Q), [n n]))
assert(isequal(size(models{3}.Q), [n n]))

R = diag(sigma_M.^2);
models{1}.R = R;
models{2}.R = R;
models{3}.R = R;
assert(isequal(size(R), [ny ny]))

P0 = 1000*eye(n);
%P0_init = repmat({P0}, 1, 3);
x0 = zeros(n,1);
y0 = models{1}.C * x0;
r0 = [1 1]';  % initial system mode

% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
seq = repmat({ones(1, nT+1)}, 4, 1);
seq{2}(t == t_shock(1)) = 2;  % shock 1
seq{3}(t == t_shock(2)) = 3;  % shock 2
seq{4}(t == t_shock(1)) = 2;  % both
seq{4}(t == t_shock(2)) = 3;  % both

% Define switching Kalman filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% the test simulation (t = t_shock)
SKF = SKFObserverS(models,P0,seq{4},"SKF");

% Define custom MKF test observers

% Multiple model observer 1 - with 4 sequences
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Build probability transition matrix
p_rk_g_rkm1 = [1-epsilon epsilon]';
Z = [1 1; 2 1; 1 2];  % combinations
p_rk_g_rkm1 = prod(prob_rk(Z', p_rk_g_rkm1), 1)';
p_rk_g_rkm1 = p_rk_g_rkm1 ./ sum(p_rk_g_rkm1);  % normalized
T = repmat(p_rk_g_rkm1', 3, 1);

MKF1 = MKFObserverS(models,P0,seq,T,"MKF1");
assert(MKF1.nh == size(seq, 1))

% Multiple model observer 2 - has only the correct sequence
seq = {ones(1, nT+1)};
seq{1}(t == t_shock(1)) = 2;
seq{1}(t == t_shock(2)) = 3;
MKF2 = MKFObserverS(models,P0,seq,T,"MKF2");
assert(MKF2.nh == 1)

% Multiple model observer with sequence fusion
SF1_seq = { ...
    [1 1 1]
    [2 1 1]
    [1 2 1]
    [1 1 2]
    [3 1 1]
    [1 3 1]
    [1 1 3] ...
};
MKF_SF = MKFObserverSF(models,P0,SF1_seq,T,"MKF_SF");
assert(MKF_SF.nh == 13)
assert(isequal(MKF_SF.r0, [1 1 1 1 1 1 1 2 1 1 1 1 3]'))

% Multiple model observer with sequence fusion and detection
% interval > 1
d = 3;
MKF_SF_DI = MKFObserverSF_DI(models,P0,SF1_seq,T,d,"MKF_SF_DI");
assert(MKF_SF_DI.nh == 13)
assert(isequal(MKF_SF_DI.r0, [1 1 1 1 1 1 1 2 1 1 1 1 3]'))
assert(MKF_SF_DI.id == 0)
assert(MKF_SF_DI.id_next == 1)

% Multiple model observer with sequence pruning
nh = 9;
n_min = 3;
MKF_SP = MKFObserverSP(models,P0,T,nh,n_min,"MKF_SP");
assert(MKF2.nh == 1)

% Multiple model observer with sequence pruning and detection
% interval > 1
d = 3;
nh = 9;
n_min = 3;
MKF_SP_DI = MKFObserverSP_DI(models,P0,T,d,nh,n_min,'MKF_SP_DI');
assert(MKF2.nh == 1)

% Choose observers to test
observers = {SKF, MKF1, MKF2, MKF_SF, MKF_SF_DI, MKF_SP, MKF_SP_DI};
%observers = {SKF, MKF_SF, MKF_SF_DI};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = find(obs_labels == "MKF_SF_DI");

% Simulate observers
[Xk_est,Yk_est,DiagP,MKF_vars] = ...
    run_simulation_obs(Ym,U,alpha,seq,observers,f_mkf);

% Output estimation errors
E_obs = repmat(Y, 1, n_obs) - Yk_est;

% Combine and display results
sim_results = table(t,alpha,U,X,Y,Ym,Xk_est,Yk_est,E_obs);
% disp(sim_results)

% % Plot observer estimates
% figure(5); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,escape_latex_chars(obs_labels))

% Plot MKF observer variables
% figure(6); clf
% 
% if isfield(MKF_vars, "merged")
%     TrP = MKF_vars.merged.trP_obs;
%     P_seq_g_Yk = cell2mat(MKF_vars.merged.p_seq_g_Yk);
%     MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.merged.Yk_est, ...
%         'UniformOutput', false));
% else
%     TrP = MKF_vars.trP_obs;
%     P_seq_g_Yk = MKF_vars.p_seq_g_Yk;
%     MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.Yk_est, ...
%         'UniformOutput', false));
% end
% nh_plot = size(P_seq_g_Yk, 2);
% 
% ax1 = subplot(2,1,1);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Tr}(\mathbf{P}(k))$'};
% make_waterfall_plot(t, TrP, [0 round(max(TrP(10:end,:),[],[1 2]))], ...
%     ax_labels, [0 64]);
% title("Trace of error covariance")
% 
% ax2 = subplot(2,1,2);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Pr}(\Gamma_f(k) | \mathbf{Y}(k))$'};
% make_waterfall_plot(t, P_seq_g_Yk, [0 1], ax_labels, [0 64]);
% title("Hypothesis probabilities")
% 
% linkprop([ax1 ax2], 'xlim');


% Plot MKF observer variables
% figure(7); clf
% 
% ax1 = subplot(2,1,1);
% plot(t, Y(:, 1), '--'); hold on
% plot(t, MKF_Yk_est(:, 1:2:nh_plot*2));
% plot(t, Yk_est(:, 1), 'k.-');
% grid on
% ylabel("$y_1(k)$", 'Interpreter', 'latex')
% legend(["$y_1(k)$" compose("$%s{y}_{%d,1}(k)$", "\hat", 1:nh_plot) ...
%     sprintf("$%s{y_1}(k)$", "\hat")], 'Interpreter', 'latex')
% 
% ax2 = subplot(2,1,2);
% plot(t, Y(:, 2), '--'); hold on
% plot(t, MKF_Yk_est(:, (1:2:nh_plot*2) + 1));
% plot(t, Yk_est(:, 2), 'k.-');
% grid on
% ylabel("$y_2(k)$", 'Interpreter', 'latex')
% legend(["$y_2(k)$" compose("$%s{y}_{%d,2}(k)$", "\hat", 1:nh_plot) ...
%     sprintf("$%s{y_2}(k)$", "\hat")], 'Interpreter', 'latex')
% xlabel('Time ($t$)', 'Interpreter', 'latex')
% 
% linkaxes([ax1 ax2], 'x')

% Compute mean-squared errors
MSE = struct();
for i = 1:n_obs
    MSE.(observers{i}.label) = mean(E_obs(:, i*ny-1:i*ny).^2);
end
% disp(MSE)

% Check final states and outputs of all observers
XYk_est = cell2table(cellfun(@(obs) [obs.xk_est' obs.yk_est'], ...
    observers', 'UniformOutput', false), ...
    'RowNames', obs_labels);
test_XYk_est = array2table( ...
    [ ...
   -1.7921    8.8475    0.9998    0.9997   -0.1989    0.9821
   -1.7921    8.8475    0.9998    0.9997   -0.1989    0.9821
   -1.7921    8.8475    0.9998    0.9997   -0.1989    0.9821
   -1.7911    8.8493    1.0001    1.0002   -0.1988    0.9823
   -1.7924    8.8481    0.9997    0.9998   -0.1990    0.9821
   -1.7864    8.8510    1.0014    1.0009   -0.1983    0.9825
   -1.7925    8.8488    0.9997    1.0000   -0.1990    0.9822...
    ], 'RowNames', ...
    {'SKF', 'MKF1', 'MKF2', 'MKF_SF', 'MKF_SF_DI', 'MKF_SP', 'MKF_SP_DI'} ...
);
assert(isequal(round(XYk_est.Variables, 4), ...
    test_XYk_est{XYk_est.Properties.RowNames, :}))

% Check final error covariance estimates
final_DiagP = table(reshape(DiagP(t == t(end), :), [], n_obs)', ...
    'RowNames', obs_labels, 'VariableNames', {'DiagP_final'});
test_final_DiagP = table( ...
    [ ...
    0.0405    0.0407    0.0015    0.0015
    0.0405    0.0407    0.0015    0.0015
    0.0405    0.0407    0.0015    0.0015
    0.0473    0.0474    0.0128    0.0128
    0.0412    0.0414    0.0068    0.0068
    0.0414    0.0415    0.0120    0.0120
    0.0405    0.0408    0.0015    0.0015 ...
    ], 'RowNames', ...
    {'SKF', 'MKF1', 'MKF2', 'MKF_SF', 'MKF_SF_DI', 'MKF_SP', 'MKF_SP_DI'}, ...
    'VariableNames', {'DiagP_final'} ...
);
assert(isequal(round(final_DiagP.Variables, 4), ...
    test_final_DiagP{final_DiagP.Properties.RowNames, :}))

% Display trace of covariance matrix data for MKF observer filters
% disp([table(t) array2table(MKF_vars.trP_obs, 'VariableNames', ...
%     compose("Tr(P_%d)", 1:observers{f_mkf}.nh))])

% Results on 2022-11-29
MSE_test_values = struct( ...
 'SKF', [0.000009, 0.000010],  ...
 'MKF1', [0.000228, 0.000194],  ...
 'MKF2', [0.000009, 0.000010], ...
 'MKF_SF', [0.000249, 0.000207], ...
 'MKF_SF_DI', [0.000282, 0.000213], ...
 'MKF_SP', [0.000266, 0.000177], ...
 'MKF_SP_DI', [0.000914, 0.000957] ...
);

labels = fieldnames(MSE);
% for i = 1:numel(labels)
%     fprintf("%15s: %.6f, %.6f (%.6f, %.6f)\n", labels{i}, MSE.(labels{i}), ...
%         MSE_test_values.(labels{i}))
% end
for i = 1:numel(labels)
    assert(isequal(round(MSE.(labels{i}), 6), MSE_test_values.(labels{i})))
end


% Now re-run simulation with noise

% Choose measurement noise for plant
sigma_MP = [0.05; 0.05];  % See above for simulation without noise
V = sigma_MP'.*randn(nT+1, ny);

% Run simulation
[X, Y, Ym] = run_simulation_sys(sys_models,U_sim,V,alpha,nT,x0);

% Reset observers to initial states
for f = 1:n_obs
    observers{f}.reset()
end

% Simulate observers
[Xk_est,Yk_est,DiagP,MKF_vars] = ...
    run_simulation_obs(Ym,U,alpha,seq,observers,f_mkf);

% % Plot of inputs and outputs
% figure(8); clf
% P = cumsum(Wp);
% u_labels = [compose("$u_%d(k)$", 1:size(U,2)) compose("$p_%d(k)$", 1:size(P,2))];
% y_labels = [compose("$y_%d(k)$", 1:ny) compose("$y_{m,%d}(k)$", 1:ny)];
% x_label = "Time, $t$";
% make_iodplot(Y,Ym,t,[U P],u_labels,y_labels,x_label)

% Output estimation errors
E_obs = repmat(Y, 1, n_obs) - Yk_est;

% Combine and display results
sim_results = table(t,alpha,U,X,Y,Ym,Xk_est,Yk_est,E_obs);
%disp(sim_results)

% % Plot observer estimates
% figure(9); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,escape_latex_chars(obs_labels))
% 
% % Plot MKF observer variables
% figure(10); clf
% 
% if isfield(MKF_vars, "merged")
%     TrP = MKF_vars.merged.trP_obs;
%     P_seq_g_Yk = cell2mat(MKF_vars.merged.p_seq_g_Yk);
%     MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.merged.Yk_est, ...
%         'UniformOutput', false));
% else
%     TrP = MKF_vars.trP_obs;
%     P_seq_g_Yk = MKF_vars.p_seq_g_Yk;
%     MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.Yk_est, ...
%         'UniformOutput', false));
% end
% nh_plot = size(P_seq_g_Yk, 2);
% 
% ax1 = subplot(2,1,1);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Tr}(\mathbf{P}(k))$'};
% make_waterfall_plot(t, TrP, [0 ceil(max(TrP(10:end,:),[],[1 2]))], ...
%     ax_labels, [0 64]);
% title("Trace of error covariance")
% 
% ax2 = subplot(2,1,2);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Pr}(\Gamma_f(k) | \mathbf{Y}(k))$'};
% make_waterfall_plot(t, P_seq_g_Yk, [0 1], ax_labels, [0 64]);
% title("Hypothesis probabilities")
% 
% linkprop([ax1 ax2], 'xlim');
% 
% MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.Yk_est, ...
%     'UniformOutput', false));
% 
% % Plot MKF observer variables
% figure(11); clf
% 
% ax1 = subplot(2,1,1);
% plot(t, Y(:, 1), '--'); hold on
% plot(t, MKF_Yk_est(:, 1:2:nh_plot*2));
% plot(t, Yk_est(:, 1), 'k.-');
% grid on
% ylabel("$y_1(k)$", 'Interpreter', 'latex')
% legend(["$y_1(k)$" compose("$%s{y}_{%d,1}(k)$", "\hat", 1:nh_plot) ...
%     sprintf("$%s{y_1}(k)$", "\hat")], 'Interpreter', 'latex')
% 
% ax2 = subplot(2,1,2);
% plot(t, Y(:, 2), '--'); hold on
% plot(t, MKF_Yk_est(:, (1:2:nh_plot*2) + 1));
% plot(t, Yk_est(:, 2), 'k.-');
% grid on
% ylabel("$y_2(k)$", 'Interpreter', 'latex')
% legend(["$y_2(k)$" compose("$%s{y}_{%d,2}(k)$", "\hat", 1:nh_plot) ...
%     sprintf("$%s{y_2}(k)$", "\hat")], 'Interpreter', 'latex')
% xlabel('Time ($t$)', 'Interpreter', 'latex')
% 
% linkaxes([ax1 ax2], 'x')


% Compute mean-squared errors
MSE = struct();
for i = 1:n_obs
    MSE.(observers{i}.label) = mean(E_obs(:, i*ny-1:i*ny).^2);
end
% disp(MSE);

% Check final states and outputs of all observers
XYk_est = cell2table(cellfun(@(obs) [obs.xk_est' obs.yk_est'], ...
    observers', 'UniformOutput', false), ...
    'RowNames', obs_labels);
test_XYk_est = array2table( ...
    [ ...
   -1.8946    8.7429    0.9832    0.9827   -0.2103    0.9705
   -1.8946    8.7429    0.9832    0.9827   -0.2103    0.9705
   -1.8946    8.7429    0.9832    0.9827   -0.2103    0.9705
   -1.8800    8.7599    0.9916    0.9916   -0.2087    0.9724
   -1.8940    8.7285    0.9826    0.9792   -0.2102    0.9689
   -1.8815    8.7508    0.9873    0.9858   -0.2089    0.9713
   -1.8935    8.7283    0.9825    0.9789   -0.2102    0.9688 ...
    ], 'RowNames', ...
    {'SKF', 'MKF1', 'MKF2', 'MKF_SF', 'MKF_SF_DI', 'MKF_SP', 'MKF_SP_DI'} ...
);
assert(isequal(round(XYk_est.Variables, 4), ...
    test_XYk_est{XYk_est.Properties.RowNames, :}))

% Check final error covariance estimates
% fprintf("%10.6f %10.6f %10.6f %10.6f\n", final_DiagP.Variables')
final_DiagP = table(reshape(DiagP(t == t(end), :), [], n_obs)', ...
    'RowNames', obs_labels, 'VariableNames', {'DiagP_final'});
test_final_DiagP = table( ...
    [ ...
  0.040535   0.040672   0.001508   0.001516
  0.040535   0.040672   0.001508   0.001516
  0.040535   0.040672   0.001508   0.001516
  0.052138   0.052132   0.013692   0.013692
  0.041292   0.041557   0.006921   0.007005
  0.041397   0.041535   0.012044   0.012063
  0.040570   0.040800   0.001512   0.001527 ...
    ], 'RowNames', ...
    {'SKF', 'MKF1', 'MKF2', 'MKF_SF', 'MKF_SF_DI', 'MKF_SP', 'MKF_SP_DI'}, ...
    'VariableNames', {'DiagP_final'} ...
);
assert(isequal(round(final_DiagP.Variables, 6), ...
    test_final_DiagP{final_DiagP.Properties.RowNames, :}))

% Display trace of covariance matrix data for MKF observer filters
% disp([table(t) array2table(MKF_vars.trP_obs, 'VariableNames', ...
%     compose("Tr(P_%d)", 1:observers{f_mkf}.nh))])

% Results on 2022-11-29
% Note: these results may depend somewhat on the random initialization
MSE_test_values = struct( ...
 'SKF', [0.000606, 0.000692],  ...
 'MKF1', [0.001051, 0.000828],  ...
 'MKF2', [0.000606, 0.000692],  ...
 'MKF_SF', [0.001182, 0.001006], ...
 'MKF_SF_DI', [0.001098, 0.001039], ...  % was [0.001703, 0.000884] before fixing prior prob.
 'MKF_SP', [0.001169, 0.000898], ...
 'MKF_SP_DI', [0.001707 0.000882] ...
);

labels = fieldnames(MSE);
% for i = 1:numel(labels)
%     fprintf("%15s: %.6f, %.6f (%.6f, %.6f)\n", labels{i}, MSE.(labels{i}), ...
%         MSE_test_values.(labels{i}))
% end
for i = 1:numel(labels)
    assert(isequal(round(MSE.(labels{i}), 6), MSE_test_values.(labels{i})))
end



%% Simulation test on alternative 2x2 linear system with RODD disturbances

% Set rng here - better for debugging
%rng(0)

% Sample time
Ts = 1;

% Discrete time state space model
A = [0.8890         0    1  0.5
          0    0.8890  0.5    1
          0         0    1    0
          0         0    0    1];
B = [    1  0.5  0  0;
       0.5    1  0  0;
         0    0  1  0;
         0    0  0  1];
C = [ 0.1110       0  0  0;
           0  0.1110  0  0];
D = zeros(2, 4);
Gpss = ss(A,B,C,D,Ts);

% Dimensions
[n, nu, ny] = check_dimensions(A, B, C, D);

% Simulation settings
nT = 50;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = [5 15];
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
%P = cumsum(Wp);
U_sim = [U Wp];
x0 = zeros(n, 1);

% System models
sys_model.A = A;
sys_model.B = B;  % input disturbances unmeasured
sys_model.C = C;
sys_model.Ts = Ts;
sys_models = repmat({sys_model}, 1, 3);

% Choose measurement noise for plant
sigma_MP = [0; 0];
V = sigma_MP'.*randn(nT+1, ny);

% Run simulation
[X, Y, Ym] = run_simulation_sys(sys_models,U_sim,V,alpha,nT,x0);

% % Simulate system
% X2 = zeros(nT+1,n);
% Y2 = zeros(nT+1,ny);
% xk = x0;
% for i = 1:nT+1
% 
%     % Inputs
%     uk = U_sim(i,:)';
% 
%     % Compute y(k)
%     yk = C * xk + D * uk;
% 
%     % Store results
%     X2(i, :) = xk';
%     Y2(i, :) = yk';
% 
%     % Compute x(k+1)
%     xk = A * xk + B * uk;
% 
% end
% 
% % Check simulation output is correct
% [Y3, t, X3] = lsim(Gpss, U_sim, t, x0);
% assert(isequal(X, X2))
% assert(isequal(Y, Y2))
% assert(isequal(X, X3))
% assert(isequal(Y, Y3))

% Plot of inputs and outputs
% figure(4); clf
% P = cumsum(Wp);
% u_labels = [compose("$u_%d(k)$", 1:size(U,2)) compose("$p_%d(k)$", 1:size(P,2))];
% y_labels = [compose("$y_%d(k)$", 1:ny) compose("$y_{m,%d}(k)$", 1:ny)];
% x_label = "Time, $t$";
% make_iodplot(Y,Ym,t,[U P],u_labels,y_labels,x_label)

% Designate measured input and output signals
u_known = [true; true; false; false];
y_meas = [true; true];

% Observer model (without unmeasured disturbance input)
Bu = B(:, u_known);
Du = D(:, u_known);
nu = sum(u_known);
nw = sum(~u_known);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_known);
nw = sum(~u_known);

% Dimensions of observer model
[n, nu, ny] = check_dimensions(A, Bu, C, Du);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.05; 0.05];
%sigma_M = [0; 0];  % set to zero for testing
sigma_wp = [0.01 1; 0.01 1];

% Observer models
model.A = A;
model.B = Bu;  % input disturbances unmeasured
model.C = C;
model.Ts = Ts;
models = repmat({model}, 1, 3);

% Observer parameters (same for all observers)
models{1}.Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
models{2}.Q = diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]);
models{3}.Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2]);
assert(isequal(size(models{1}.Q), [n n]))
assert(isequal(size(models{2}.Q), [n n]))
assert(isequal(size(models{3}.Q), [n n]))

R = diag(sigma_M.^2);
models{1}.R = R;
models{2}.R = R;
models{3}.R = R;
assert(isequal(size(R), [ny ny]))

P0 = 1000*eye(n);
%P0_init = repmat({P0}, 1, 3);
x0 = zeros(n,1);
y0 = models{1}.C * x0;
r0 = [1 1]';  % initial system mode

% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
seq = repmat({ones(1, nT+1)}, 4, 1);
seq{2}(t == t_shock(1)) = 2;  % shock 1
seq{3}(t == t_shock(2)) = 3;  % shock 2
seq{4}(t == t_shock(1)) = 2;  % both
seq{4}(t == t_shock(2)) = 3;  % both

% Define switching Kalman filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% the test simulation (t = t_shock)
SKF = SKFObserverS(models,P0,seq{4},"SKF");

% Define custom MKF test observers

% Multiple model observer 1 - with 4 sequences
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Build probability transition matrix
p_rk_g_rkm1 = [1-epsilon epsilon]';
Z = [1 1; 2 1; 1 2];  % combinations
p_rk_g_rkm1 = prod(prob_rk(Z', p_rk_g_rkm1), 1)';
p_rk_g_rkm1 = p_rk_g_rkm1 ./ sum(p_rk_g_rkm1);  % normalized
T = repmat(p_rk_g_rkm1', 3, 1);

MKF1 = MKFObserverS(models,P0,seq,T,"MKF1");
assert(MKF1.nh == size(seq, 1))

% Multiple model observer 2 - has only the correct sequence
seq = {ones(1, nT+1)};
seq{1}(t == t_shock(1)) = 2;
seq{1}(t == t_shock(2)) = 3;
MKF2 = MKFObserverS(models,P0,seq,T,"MKF2");
assert(MKF2.nh == 1)

% Multiple model observer with sequence fusion
SF1_seq = { ...
    [1 1 1]
    [2 1 1]
    [1 2 1]
    [1 1 2]
    [3 1 1]
    [1 3 1]
    [1 1 3] ...
};
MKF_SF = MKFObserverSF(models,P0,SF1_seq,T,"MKF_SF");
assert(MKF_SF.nh == 13)
assert(isequal(MKF_SF.r0, [1 1 1 1 1 1 1 2 1 1 1 1 3]'))

% Multiple model observer with sequence fusion and detection
% interval > 1
d = 3;
MKF_SF_DI = MKFObserverSF_DI(models,P0,SF1_seq,T,d,"MKF_SF_DI");
assert(MKF_SF_DI.nh == 13)
assert(isequal(MKF_SF_DI.r0, [1 1 1 1 1 1 1 2 1 1 1 1 3]'))
assert(MKF_SF_DI.id == 0)
assert(MKF_SF_DI.id_next == 1)

% Multiple model observer with sequence pruning
nh = 9;
n_min = 3;
MKF_SP = MKFObserverSP(models,P0,T,nh,n_min,"MKF_SP");
assert(MKF2.nh == 1)

% Multiple model observer with sequence pruning and detection
% interval > 1
d = 3;
nh = 9;
n_min = 3;
MKF_SP_DI = MKFObserverSP_DI(models,P0,T,d,nh,n_min,'MKF_SP_DI');
assert(MKF2.nh == 1)

% Choose observers to test
observers = {SKF, MKF1, MKF2, MKF_SF, MKF_SF_DI, MKF_SP, MKF_SP_DI};
%observers = {SKF, MKF_SF, MKF_SF_DI};
n_obs = numel(observers);
obs_labels = cellfun(@(x) x.label, observers, 'UniformOutput', true);

% Identify which observer to log MKF data for
f_mkf = find(obs_labels == "MKF_SF_DI");

% Simulate observers
[Xk_est,Yk_est,DiagP,MKF_vars] = ...
    run_simulation_obs(Ym,U,alpha,seq,observers,f_mkf);

% Output estimation errors
E_obs = repmat(Y, 1, n_obs) - Yk_est;

% Combine and display results
%sim_results = table(t,alpha,U,X,Y,Ym,Xk_est,Yk_est,E_obs);
% disp(sim_results)

% % Plot observer estimates
% figure(5); clf
% plot_obs_estimates(t,X,Xk_est,Y,Yk_est,escape_latex_chars(obs_labels))

% Plot MKF observer variables
% figure(6); clf
% 
% if isfield(MKF_vars, "merged")
%     TrP = MKF_vars.merged.trP_obs;
%     P_seq_g_Yk = cell2mat(MKF_vars.merged.p_seq_g_Yk);
%     MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.merged.Yk_est, ...
%         'UniformOutput', false));
% else
%     TrP = MKF_vars.trP_obs;
%     P_seq_g_Yk = MKF_vars.p_seq_g_Yk;
%     MKF_Yk_est = cell2mat(cellfun(@(x) x', MKF_vars.Yk_est, ...
%         'UniformOutput', false));
% end
% nh_plot = size(P_seq_g_Yk, 2);
% 
% ax1 = subplot(2,1,1);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Tr}(\mathbf{P}(k))$'};
% make_waterfall_plot(t, TrP, [0 round(max(TrP(10:end,:),[],[1 2]))], ...
%     ax_labels, [0 64]);
% title("Trace of error covariance")
% 
% ax2 = subplot(2,1,2);
% ax_labels = {'Time ($t$)', 'Hypotheses, $f$', ...
%     '$\mathrm{Pr}(\Gamma_f(k) | \mathbf{Y}(k))$'};
% make_waterfall_plot(t, P_seq_g_Yk, [0 1], ax_labels, [0 64]);
% title("Hypothesis probabilities")
% 
% linkprop([ax1 ax2], 'xlim');


% Plot MKF observer variables
% figure(7); clf
% 
% ax1 = subplot(2,1,1);
% plot(t, Y(:, 1), '--'); hold on
% plot(t, MKF_Yk_est(:, 1:2:nh_plot*2));
% plot(t, Yk_est(:, 1), 'k.-');
% grid on
% ylabel("$y_1(k)$", 'Interpreter', 'latex')
% legend(["$y_1(k)$" compose("$%s{y}_{%d,1}(k)$", "\hat", 1:nh_plot) ...
%     sprintf("$%s{y_1}(k)$", "\hat")], 'Interpreter', 'latex')
% 
% ax2 = subplot(2,1,2);
% plot(t, Y(:, 2), '--'); hold on
% plot(t, MKF_Yk_est(:, (1:2:nh_plot*2) + 1));
% plot(t, Yk_est(:, 2), 'k.-');
% grid on
% ylabel("$y_2(k)$", 'Interpreter', 'latex')
% legend(["$y_2(k)$" compose("$%s{y}_{%d,2}(k)$", "\hat", 1:nh_plot) ...
%     sprintf("$%s{y_2}(k)$", "\hat")], 'Interpreter', 'latex')
% xlabel('Time ($t$)', 'Interpreter', 'latex')
% 
% linkaxes([ax1 ax2], 'x')

% Compute mean-squared errors
MSE = struct();
for i = 1:n_obs
    MSE.(observers{i}.label) = mean(E_obs(:, i*ny-1:i*ny).^2);
end
% disp(MSE)

% Check final states and outputs of all observers
XYk_est = cell2table(cellfun(@(obs) [obs.xk_est' obs.yk_est'], ...
    observers', 'UniformOutput', false), ...
    'RowNames', obs_labels);
test_XYk_est = array2table( ...
    [ ...
    4.3994    8.8325    1.0000    0.9999    0.4883    0.9804
    4.3994    8.8325    1.0000    0.9999    0.4883    0.9804
    4.3994    8.8325    1.0000    0.9999    0.4883    0.9804
    4.4001    8.8322    1.0003    0.9997    0.4884    0.9804
    4.3990    8.8331    0.9998    1.0002    0.4883    0.9805
    4.4079    8.8240    1.0045    0.9955    0.4893    0.9795
    4.3986    8.8338    0.9995    1.0005    0.4882    0.9806...
    ], 'RowNames', ...
    {'SKF', 'MKF1', 'MKF2', 'MKF_SF', 'MKF_SF_DI', 'MKF_SP', 'MKF_SP_DI'} ...
);
assert(isequal(round(XYk_est.Variables, 4), ...
    test_XYk_est{XYk_est.Properties.RowNames, :}))

% Check final error covariance estimates
final_DiagP = table(reshape(DiagP(t == t(end), :), [], n_obs)', ...
    'RowNames', obs_labels, 'VariableNames', {'DiagP_final'});
test_final_DiagP = table( ...
    [ ...
    0.0404    0.0405    0.0019    0.0019
    0.0404    0.0405    0.0019    0.0019
    0.0404    0.0405    0.0019    0.0019
    0.0460    0.0460    0.0131    0.0131
    0.0410    0.0412    0.0068    0.0068
    0.0412    0.0413    0.0123    0.0123
    0.0404    0.0406    0.0019    0.0019 ...
    ], 'RowNames', ...
    {'SKF', 'MKF1', 'MKF2', 'MKF_SF', 'MKF_SF_DI', 'MKF_SP', 'MKF_SP_DI'}, ...
    'VariableNames', {'DiagP_final'} ...
);
assert(isequal(round(final_DiagP.Variables, 4), ...
    test_final_DiagP{final_DiagP.Properties.RowNames, :}))

% Display trace of covariance matrix data for MKF observer filters
% disp([table(t) array2table(MKF_vars.trP_obs, 'VariableNames', ...
%     compose("Tr(P_%d)", 1:observers{f_mkf}.nh))])

% Results on 2022-11-29
MSE_test_values = struct( ...
 'SKF', [0.000008, 0.000008],  ...
 'MKF1', [0.000196, 0.000185],  ...
 'MKF2', [0.000008, 0.000008], ...
 'MKF_SF', [0.000209, 0.000191], ...
 'MKF_SF_DI', [0.000269, 0.000205], ...
 'MKF_SP', [0.000198, 0.000217], ...
 'MKF_SP_DI', [0.000787, 0.000760] ...
);

labels = fieldnames(MSE);
% for i = 1:numel(labels)
%     fprintf("%15s: %.6f, %.6f (%.6f, %.6f)\n", labels{i}, MSE.(labels{i}), ...
%         MSE_test_values.(labels{i}))
% end
for i = 1:numel(labels)
    assert(isequal(round(MSE.(labels{i}), 6), MSE_test_values.(labels{i})))
end
