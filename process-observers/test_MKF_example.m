% Test simulation of Kalman filters and multi-model
% observers to estimate randomly-occurring deterministic 
% disturbances (RODDs).
%

clear all

%show_plots = true;
show_plots = false;

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

if show_plots
    figure(1)
    subplot(3,1,[1 2])
    stairs(0:nT,Wp); grid on
    ylabel('wp(k)')
    title('Disturbance')
    subplot(3,1,3)
    stairs(0:nT, alpha); grid on
    ylim([-0.1 1.1])
    xlabel('k')
    ylabel('alpha(k)')
end

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

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = KFPSS;
for i = 1:nT+1
    Xk_est(i,:) = obs.xkp1_est';
    Yk_est(i,:) = obs.ykp1_est';
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs.update(yk, uk);
end

if show_plots
    figure(3)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)','x2_est(k)')
    title('States')
end

% Calculate mean-squared error in state estimates
mse_KFPSS = mean((X - Xk_est).^2, [1 2]);

% Store results
Xk_est_KFPSS = Xk_est;
Yk_est_KFPSS = Yk_est;

% Filtering form
KFFSS = KalmanFilterFSS(obs_model,'KFFSS');

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = KFFSS;
for i = 1:nT+1
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs.update(yk, uk);
    Xk_est(i,:) = obs.xk_est';
    Yk_est(i,:) = obs.yk_est';
end

if show_plots
    figure(4)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)','x2_est(k)')
    title('States')
end

% Calculate mean-squared error in state estimates
mse_KFFSS = mean((X - Xk_est).^2, [1 2]);

% Store results
Xk_est_KFFSS = Xk_est;
Yk_est_KFFSS = Yk_est;

% Dynamic Kalman filters
P0 = eye(n);

% Prediction form
KFP = KalmanFilterP(obs_model,P0,'KFP');

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = KFP;
for i = 1:nT+1
    Xk_est(i,:) = obs.xkp1_est';
    Yk_est(i,:) = obs.ykp1_est';
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs.update(yk, uk);
end

if show_plots
    figure(5)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)','x2_est(k)')
    title('States')
end

% Calculate mean-squared error in state estimates
mse_KFP = mean((X - Xk_est).^2, [1 2]);

% Store results
Xk_est_KFP = Xk_est;
Yk_est_KFP = Yk_est;

% Filtering form
KFF = KalmanFilterF(obs_model,P0,'KFF');

% Simulate observers
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = KFF;
for i = 1:nT+1
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs.update(yk, uk);
    Xk_est(i,:) = obs.xk_est';
    Yk_est(i,:) = obs.yk_est';
end

if show_plots
    figure(6)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)','x2_est(k)')
    title('States')
end

% Calculate mean-squared error in state estimates
mse_KFF = mean((X - Xk_est).^2, [1 2]);

% Store results
Xk_est_KFF = Xk_est;
Yk_est_KFF = Yk_est;

% Switching Kalman filter - scheduled
obs_models = {obs_model, obs_model};
P0 = eye(n);
obs_models{1}.Q = diag([0.01^2 sigma_wp{1}(1)^2]);
obs_models{2}.Q = diag([0.01^2 sigma_wp{1}(2)^2]);
seq = alpha' + 1;
SKF = SKFObserverS(obs_models,P0,seq,"SKF");

% Simulate observers
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = SKF;
for i = 1:nT+1
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs.update(yk, uk);
    Xk_est(i,:) = obs.xk_est';
    Yk_est(i,:) = obs.yk_est';
end

if show_plots
    figure(7)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)','x2_est(k)')
    title('States')
end

% Calculate mean-squared error in state estimates
mse_SKF = mean((X - Xk_est).^2, [1 2]);

% Store results
Xk_est_SKF = Xk_est;
Yk_est_SKF = Yk_est;

% Sub-optimal multi-model observer 1 simulation
% Sub-optimal multi-model observer as described by Robertson _et al._ (1995).

% Multi-model observer 1 - sequence fusion
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.D = D;
model.Ts = Ts;
P0 = eye(n);
Q0 = diag([0.01^2 0]);
R = 0.1^2;
nf = 5;  % fusion horizon
m = 1;  % maximum number of shocks
d = 3;  % spacing parameter
io.u_known = u_known;
io.y_meas = y_meas;
MKF1 = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nf,m,d,'MKF1');

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = MKF1;
for i = 1:nT+1
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs.update(yk, uk);
    Xk_est(i,:) = obs.xk_est';
    Yk_est(i,:) = obs.yk_est';
end

if show_plots
    figure(8)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)',['x2_est(k)'])
    title('States')
end

% Calculate mean-squared error in state estimates
mse_MKF1 = mean((X - Xk_est).^2, [1 2]);

% Store results
Xk_est_MKF1 = Xk_est;
Yk_est_MKF1 = Yk_est;

% Sub-optimal multi-model observer 2 simulation
% Sub-optimal multi-model observer as described by Eriksson and Isaksson (1996).

% Multi-model observer 2 - sequence pruning
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.D = D;
model.Ts = Ts;
P0 = eye(n);
Q0 = diag([0.01^2 0]);
R = 0.1^2;
nh = 5;  % number of filters
n_min = 2;  % minimum life of cloned filters
io.u_known = u_known;
io.y_meas = true(ny, 1);
MKF2 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,'MKF2');

% Simulate observer
Xk_est = nan(nT+1,n);
Yk_est = nan(nT+1,ny);
obs = MKF2;
for i = 1:nT+1
    uk = U(i,:)';
    yk = Ym(i,:)';
    obs.update(yk, uk);
    Xk_est(i,:) = obs.xk_est';
    Yk_est(i,:) = obs.yk_est';
end

if show_plots
    figure(9)
    subplot(2,1,1)
    plot(t,[Y Yk_est]); grid on
    ylabel('y(k) and y_est(k)')
    legend('y(k)','y_est(k)')
    title('Output')
    subplot(2,1,2)
    plot(t, [X Xk_est]); grid on
    xlabel('k')
    ylabel('xi(k) and xi_est(k)')
    legend('x1(k)','x2(k)','x1_est(k)',['x2_est(k)'])
    title('States')
end

% Calculate mean-squared error in state estimates
mse_MKF2 = mean((X - Xk_est).^2, [1 2]);

% Store results
Xk_est_MKF2 = Xk_est;
Yk_est_MKF2 = Yk_est;

% Check MSE results for all observer
% Note: 2022-07-04 Updated result after fixing MKFObserverSF sigma_wp
assert(round(mse_KFPSS, 6) == 0.086430)
assert(round(mse_KFFSS, 6) == 0.049202)
assert(round(mse_KFP, 6) == 0.090836)
assert(round(mse_KFF, 6) == 0.052274)
assert(round(mse_SKF, 6) == 0.022272)
%assert(round(mse_MKF1, 4) == 0.1280)  % TODO: Not working yet
assert(round(mse_MKF2, 6) == 0.029665)

% Combine and save all results to csv file
% These are used by test_MKF_example_sim.m
sim_results = table(t, ...
    Xk_est_KFPSS, Yk_est_KFPSS, ...
    Xk_est_KFFSS, Yk_est_KFFSS, ...
    Xk_est_KFP, Yk_est_KFP, ...
    Xk_est_KFF, Yk_est_KFF, ...
    Xk_est_SKF, Xk_est_SKF, ...
    Xk_est_MKF1, Yk_est_MKF1, ...
    Xk_est_MKF2, Yk_est_MKF2);
filename = 'MKF_example.csv';
results_dir = 'results';
writetable(sim_results, fullfile(results_dir, filename));
