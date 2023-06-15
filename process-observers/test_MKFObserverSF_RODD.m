% Test MKFObserverSF_RODD class
%
% TODO: I think this can be deleted as there are enough tests 
% in the test_MKFObserverSF script.
% Keep the sequence generation tests though if they aren't
% somewhere else.
% Maybe copy the initialization tests from here.
% Also: tests for MKFObserverSF_RODD95 need to be moved to
% test_MKFObserverSF.m
%

clear all

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
f = 3;  % number of detection intervals in fusion horizon
m = 1;  % maximum number of shocks
d = 3;  % spacing parameter
io = struct('u_known', u_known, 'y_meas', true(ny, 1));
MKF_SF = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,"MKF_SF");

seq_test = {
    [1 1 1]; ...
    [2 1 1]; ...
    [1 2 1]; ...
    [1 1 2]
};
assert(isequal(MKF_SF.seq, seq_test))

% Check hypothesis probabilities
alpha = 1 - [1-epsilon].^d;
assert(MKF_SF.alpha == alpha)
assert(isequal(MKF_SF.p_rk, [(1-alpha) alpha]'))
p_seq = [(1-alpha).^3 alpha*(1-alpha).^2 alpha*(1-alpha).^2 alpha*(1-alpha).^2]';
assert(isequal(MKF_SF.p_seq, p_seq))
assert(MKF_SF.beta == sum(p_seq));


%% Test sequence generation 2

% Load SISO system and disturbance model from file
sys_rodin_step

% Define sequence fusion observer
P0 = eye(n);
Q0 = diag([0.01 0]);
R = sigma_M^2;
f = 2;  % number of detection intervals in fusion horizon
m = 2;  % maximum number of shocks
d = 5;  % spacing parameter
io = struct('u_known', u_known, 'y_meas', true(ny, 1));
MKF_SF = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,"MKF_SF");

seq_test = { ...
    [1 1]; ...
    [2 1]; ...
    [1 2]; ...
    [2 2] ...
};
assert(isequal(MKF_SF.seq, seq_test))

% Check hypothesis probabilities
alpha = 1 - [1-epsilon].^d;
assert(MKF_SF.alpha == alpha)
assert(isequal(MKF_SF.p_rk, [(1-alpha) alpha]'))
p_seq = [(1-alpha).^2 alpha*(1-alpha) alpha*(1-alpha) alpha.^2]';
assert(max(abs(MKF_SF.p_seq - p_seq)) < 1e-15)
assert(MKF_SF.beta == sum(p_seq));


%% Test observer initialization - SISO system

% Load SISO system and disturbance model from file
sys_rodin_step
u_known = u_known;  % variable name changed

% Load observers from file
obs_rodin_step

test_seq = { ...
    [1 1 1]; ...
    [2 1 1]; ...
    [1 2 1]; ...
    [1 1 2]
};

% Check observer attributes
%MKF_SF1 = MKFObserverSF_RODD(model,io,P0,epsilon, ...
%    sigma_wp{1},Q0,R,nf,m,d,label);
assert(strcmp(MKF_SF1.type, "MKF_SF_RODD"))
assert(isequal(MKF_SF1.sys_model, model))
assert(isequaln(MKF_SF1.io, io))
assert(isequal(MKF_SF1.P0, P0))
assert(MKF_SF1.epsilon == epsilon)
assert(isequal(MKF_SF1.sigma_wp, sigma_wp))
assert(isequal(MKF_SF1.Q0, Q0))
assert(isequal(MKF_SF1.R, R))
assert(MKF_SF1.nf == 3)
assert(MKF_SF1.m == 1)
assert(MKF_SF1.d == 5)

assert(MKF_SF1.n == 2)
assert(MKF_SF1.nu == 1)
assert(MKF_SF1.ny == 1)
assert(MKF_SF1.nj == 2)
assert(MKF_SF1.nh == 6)
assert(MKF_SF1.nm == 4)
assert(MKF_SF1.nh_max == 6)
assert(isequaln([MKF_SF1.i MKF_SF1.i_next], [3 1]))
assert(isequaln([MKF_SF1.id MKF_SF1.id_next], [0 1]))
assert(MKF_SF1.Ts == model.Ts)
assert(isequal(MKF_SF1.models{1}.Q, [0.01 0; 0 sigma_wp{1}(1)^2]))
assert(isequal(MKF_SF1.models{2}.Q, [0.01 0; 0 sigma_wp{1}(2)^2/MKF_SF1.d]))
assert(isequal(MKF_SF1.R, R))
assert(isequal(MKF_SF1.models{1}.R, R))
assert(isequal(MKF_SF1.models{2}.R, R))
assert(isequal(size(MKF_SF1.filters.Xkp1_est), [MKF_SF1.n 1 MKF_SF1.nh]))
assert(isequal(size(MKF_SF1.filters.Pkp1), [MKF_SF1.n MKF_SF1.n MKF_SF1.nh]))
assert(isequal(size(MKF_SF1.seq), [MKF_SF1.nm 1]))
assert(isequal(size(cell2mat(MKF_SF1.seq)), [MKF_SF1.nm MKF_SF1.nf]))
assert(isequal(MKF_SF1.seq, test_seq))
assert(MKF_SF1.beta == sum(MKF_SF1.p_seq))
assert(isequal(MKF_SF1.xkp1_est, zeros(n,1)))
assert(isequal(MKF_SF1.Pkp1, P0))
alpha = (1 - (1 - MKF_SF1.epsilon).^MKF_SF1.d);  % prob. over detection interval 
p_rk = [1-alpha'; alpha'];
assert(isequal(round(alpha, 4), 0.0490))
assert(isequal(round(p_rk, 4), [0.9510; 0.0490]))
assert(isequal(round(MKF_SF1.alpha, 4), round(alpha, 4)))
assert(isequal(MKF_SF1.rk, [1 1 1 1 1 2]'))
assert(isequaln(MKF_SF1.p_rk_g_rkm1, nan(MKF_SF1.nh, 1)))
assert(isequaln(MKF_SF1.xk_est, nan(2, 1)))
assert(isequaln(MKF_SF1.Pk, nan(2, 2)))
assert(isequaln(MKF_SF1.yk_est, nan))

test_seq = { ...
    [1 1 1 1 1];
    [2 1 1 1 1];
    [1 2 1 1 1];
    [1 1 2 1 1];
    [1 1 1 2 1];
    [1 1 1 1 2];
    [2 2 1 1 1];
    [2 1 2 1 1];
    [2 1 1 2 1];
    [2 1 1 1 2];
    [1 2 2 1 1];
    [1 2 1 2 1];
    [1 2 1 1 2];
    [1 1 2 2 1];
    [1 1 2 1 2];
    [1 1 1 2 2]
};

% Check observer attributes
assert(strcmp(MKF_SF2.type, "MKF_SF_RODD"))
assert(isequal(MKF_SF2.sys_model, model))
assert(isequaln(MKF_SF2.io, io))
assert(isequal(MKF_SF2.P0, P0))
assert(MKF_SF2.epsilon == epsilon)
assert(isequal(MKF_SF2.sigma_wp, sigma_wp))
assert(isequal(MKF_SF2.Q0, Q0))
assert(isequal(MKF_SF2.R, R))
assert(MKF_SF2.nf == 5)
assert(MKF_SF2.m == 2)
assert(MKF_SF2.d == 3)
assert(MKF_SF2.n == 2)
assert(MKF_SF2.nu == 1)
assert(MKF_SF2.ny == 1)
assert(MKF_SF2.nj == 2)
assert(MKF_SF2.nh == 26)
assert(MKF_SF2.nm == 16)
assert(MKF_SF2.nh_max == 26)
assert(isequaln([MKF_SF2.i MKF_SF2.i_next], [5 1]))
assert(isequaln([MKF_SF2.id MKF_SF2.id_next], [0 1]))
assert(MKF_SF1.Ts == model.Ts)
assert(isequal(MKF_SF2.models{1}.Q, [0.01 0; 0 sigma_wp{1}(1)^2]))
assert(isequal(MKF_SF2.models{2}.Q, [0.01 0; 0 sigma_wp{1}(2)^2/MKF_SF2.d]))
assert(isequal(MKF_SF2.R, R))
assert(isequal(MKF_SF2.models{1}.R, R))
assert(isequal(MKF_SF2.models{2}.R, R))
assert(isequal(size(MKF_SF2.filters.Xkp1_est), [MKF_SF2.n 1 MKF_SF2.nh]))
assert(isequal(size(MKF_SF2.filters.Pkp1), [MKF_SF2.n MKF_SF2.n MKF_SF2.nh]))
assert(isequal(size(MKF_SF2.seq), [MKF_SF2.nm 1]))
assert(isequal(size(cell2mat(MKF_SF2.seq)), [MKF_SF2.nm MKF_SF2.nf]))
assert(isequal(MKF_SF2.seq, test_seq))
assert(MKF_SF2.beta == sum(MKF_SF2.p_seq))
assert(isequal(MKF_SF2.xkp1_est, zeros(n,1)))
assert(isequal(MKF_SF2.P0, P0))
assert(isequal(MKF_SF2.Pkp1, P0))
alpha = (1 - (1 - MKF_SF2.epsilon).^MKF_SF2.d);  % prob. over detection interval 
p_rk = [1-alpha'; alpha'];
assert(isequal(round(alpha, 4), 0.0297))
assert(isequal(round(p_rk, 4), [0.9703; 0.0297]))
assert(isequal(round(MKF_SF2.alpha, 4), round(alpha, 4)))
assert(isequal(MKF_SF2.rk, [ ...
    1 1 1 1 1 1 1 1 1 1 2 2 1 ...
    1 1 1 1 1 2 2 1 1 2 1 2 2]' ...
))
assert(isequaln(MKF_SF2.p_rk_g_rkm1, nan(MKF_SF2.nh, 1)))
assert(isequaln(MKF_SF2.xk_est, nan(2, 1)))
assert(isequaln(MKF_SF2.Pk, nan(2, 2)))
assert(isequaln(MKF_SF2.yk_est, nan))

% Check optional definition with an initial state estimate
label = 'MKF_testx0';
x0 = [0.1; 0.5];
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
R = sigma_M^2;
f = 15;  % fusion horizon
m = 2;  % maximum number of shocks
d = 3;  % spacing parameter
MKF_testx0 = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0);
assert(isequal(MKF_testx0.xkp1_est, x0))


%% Test observer on SISO system with 1 shock

% This test is a copy of that used for the test simulations in
% thesis report but without plot generation. See this script for
% the latest version:
%  - disturb-models/robertson-1995/rod_obs_test_sim.m

% Load SISO system and disturbance model from file
sys_rodin_step

% Simulation settings
nT = 100;
t = Ts*(0:nT)';

% Multiple model observer with sequence fusion
P0 = eye(n);
Q0 = diag([0.01 0]);
R = sigma_M^2;
f = 15;  % fusion horizon
m = 1;  % maximum number of shocks
d = 5;  % spacing parameter
label = 'MKF_SF98';
io = struct('u_known', u_known, 'y_meas', true(ny, 1));
MKF_SF98 = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% MKF_SF95 with same parameters as MKF_SF98
% TODO: Test this one
label = 'MKF_SF95';
MKF_SF95 = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model observer with sequence pruning
f = nT;  % sequence history length
nh = 5;  % number of filters
n_min = 2;  % minimum life of cloned filters
label = 'MKF_SP';
io.u_known = u_known;
io.y_meas = true(ny, 1);
MKF_SP = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,label);

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

% Custom MKF test observer
% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t == t_shock)

% Observer model without disturbance noise input
Bu = B(:, u_known);
Du = D(:, u_known);

% Multiple model filter - two sequences, one empty, one correct
obs_model = model;
obs_model.B = Bu;
obs_model.D = Du;
obs_model.R = sigma_M.^2;
obs_models = {obs_model, obs_model};
Q0 = diag([0.01 1]);
%P0_init = repmat({P0}, 1, 2);
obs_models{1}.Q = diag([Q0(1,1) sigma_wp{1}(1,1)^2]);
obs_models{2}.Q = diag([Q0(1,1) sigma_wp{1}(1,2)^2]);
seq = {ones(1, nT+1); ones(1, nT+1)};
seq{2}(t == t_shock) = 2;
p_rk = [1-epsilon epsilon]';
T = repmat(p_rk', 2, 1);
MKF3 = MKFObserverS(obs_models,P0,seq,T,'MKF3');

% Multiple model filter - one sequence with correct shock
seq = {ones(1, nT+1)};
seq{1}(t == t_shock) = 2;
p_rk = [1-epsilon epsilon]';
T = repmat(p_rk', 2, 1);
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
observers = {MKF3, MKF4, SKF, MKF_SF95, MKF_SF98, MKF_SP};

% Note: KF1 is too slow to pass static error test here

% Combine all input signals
U_sim = [U Wp];

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

% Simulate observers

% Choose measurement noise for plant simulation
sigma_MP = 0;  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(size(Y));

% Measured inputs (not including disturbances)
U_m = U;

n_obs = numel(observers);
sim_results = struct();
RMSE_results = struct();
show_plots = false;
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_result] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs, ...
        show_plots);

    % Check observer errors are zero prior to
    % input disturbance
    %assert(all(abs(sim_result.X_est(1:20,:) - X(1:20, :)) < 1e-10, [1 2]))
    %assert(all(abs(sim_result.Y_est(1:20,:) - Y(1:20, :)) < 1e-10))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
%     if all(sigma_MP == 0)
%         assert(abs(sim_result.Y_est(end, :) - Y(end, :)) < 1e-3);
%         assert(abs(sim_result.X_est(end, 2) - du0) < 1e-3);
%     end

    % Compute mean-squared error
    Y_est = sim_result.Y_est;
    RMSE_results.(obs.label) = mean((Y_est - Y).^2);
    %fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

    % Save results
    sim_results.(obs.label) = sim_result;

    % Save updated observer (not needed if using handle objects)
    observers{i} = obs;

end

% Earlier test results (with a shock of amplitude 1)
% MSE_test_values = containers.Map(...
%     {'KF2',   'KF3',   'MKF_SF1',  'MKF_SF2',  'MKF3',  'MKF4',  "SKF"}, ...
%     [0.000934 0.003524 0.004914 0.005016 0.002709 0.000929 0.000929] ...
% );

% Results on Nov 8 before reverting back the Bayesian updating
%MSE_test_values = containers.Map(...
%  {'KF2',   'KF3',   'MKF_SF1',  'MKF_SF2',  'MKF3',  'MKF4',  "SKF"}, ...
%  [0.000934 0.003524 0.009456 0.005016 0.002709 0.000929 0.000929] ...
%);
% Changes since previous results: 
%  - f, m, d parameters for MK1 changed.
%  - Shock probability and variance modified to reflect detection
%    intervals.
%  - Bayesian prob updates only at end of detection intervals
% Note: performance of MKF3 increases if shock amplitude is increased.

% Results on 2022-11-29 after modifying MKF_SF_DI
MSE_test_values = struct( ...
    'KF2', 0.000934, ...
    'KF3', 0.003524, ...
    'MKF_SF95', 0.000729, ...
    'MKF_SF98', 0.000671, ...  % 2022-12-07 after fixing MKFObserverSF_DI
    'MKF_SP', 0.000488, ...
    'MKF3', 0.000499, ...
    'MKF4', 0.000012, ...
    'SKF', 0.000012 ...
);

% for label = fieldnames(RMSE_results)'
%     fprintf("%s: %f (%f)\n", label{1}, RMSE_results.(label{1}), MSE_test_values.(label{1}))
% end
for label = fieldnames(RMSE_results)'
    assert(isequal(round(RMSE_results.(label{1}), 6), MSE_test_values.(label{1})))
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
MKF_SF1 = MKFObserverSF_RODD(sys_model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Multiple model filter 2
label = 'MKF_SF2';
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
f = 10;  % fusion horizon
m = 2;  % maximum number of shocks
d = 2;  % spacing parameter
MKF_SF2 = MKFObserverSF_RODD(sys_model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label);

% Check observer initialization
assert(isequal(MKF_SF1.epsilon, epsilon))
assert(isequal(MKF_SF1.sigma_wp, sigma_wp))
assert(MKF_SF1.nh == 19)
assert(isequaln([MKF_SF1.i MKF_SF1.i_next], [6 1]))
assert(isequaln([MKF_SF1.id MKF_SF1.id_next], [0 1]))
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
assert(isequaln(MKF_SF1.io, io))
assert(isequal(MKF_SF1.Q0, Q0))
assert(isequal(MKF_SF1.models{1}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(1).^2 sigma_wp{2}(1).^2])))
assert(isequal(MKF_SF1.models{2}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(1).^2 sigma_wp{2}(2)^2/MKF_SF1.d])))
assert(isequal(MKF_SF1.models{3}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(2)^2/MKF_SF1.d sigma_wp{2}(1).^2])))
assert(isequal(MKF_SF1.models{1}.R, R) && isequal(MKF_SF1.models{2}.R, R) ...
    && isequal(MKF_SF1.models{3}.R, R))
assert(isequal(size(MKF_SF1.seq), [MKF_SF1.nm 1]))
assert(isequal(size(cell2mat(MKF_SF1.seq)), [MKF_SF1.nm MKF_SF1.nf]))
assert(MKF_SF1.beta == sum(MKF_SF1.p_seq))
assert(isequal(MKF_SF1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SF1.Pkp1, 1000*eye(4)))
assert(sum(MKF_SF1.p_rk) == 1)
alpha = (1 - (1 - epsilon).^d);  % prob. over detection interval 
p_rk = [1-alpha'; alpha'];
Z = [0 0; 0 1; 1 0];  % combinations
p_rk = prod(prob_gamma(Z', p_rk), 1)';
p_rk = p_rk ./ sum(p_rk);  % normalized
assert(isequal(round(p_rk, 6), [0.960977; 0.019512; 0.019512]))
assert(isequal(round(MKF_SF1.p_rk, 6), [0.960977; 0.019512; 0.019512]))

% Check observer initialization
assert(isequal(MKF_SF2.epsilon, epsilon))
assert(isequal(MKF_SF2.sigma_wp, sigma_wp))
assert(MKF_SF2.nh == 331)
assert(isequaln([MKF_SF2.i MKF_SF2.i_next], [10 1]))
assert(isequaln([MKF_SF2.id MKF_SF2.id_next], [0 1]))
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
assert(isequaln(MKF_SF2.io, io))
assert(isequal(MKF_SF2.Q0, Q0))
assert(isequal(MKF_SF2.models{1}.Q, diag([0.01 0.01 0.0001 0.0001])))
assert(isequal(MKF_SF2.models{2}.Q, diag([0.01 0.01 0.0001 0.5000])))
assert(isequal(MKF_SF2.models{3}.Q, diag([0.01 0.01 0.5000 0.0001])))
assert(isequal(MKF_SF2.models{1}.R, R) && ...
    isequal(MKF_SF2.models{2}.R, R) && ...
    isequal(MKF_SF2.models{3}.R, R))
assert(isequal(size(MKF_SF2.seq), [MKF_SF2.nm 1]))
assert(isequal(size(cell2mat(MKF_SF2.seq)), [MKF_SF2.nm MKF_SF2.nf]))
assert(MKF_SF2.beta == sum(MKF_SF2.p_seq))
assert(isequal(MKF_SF2.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SF2.Pkp1, 1000*eye(4)))
assert(sum(MKF_SF2.p_rk) == 1)
alpha = (1 - (1 - epsilon).^d);  % prob. over detection interval 
p_rk = [1-alpha'; alpha'];
Z = [0 0; 0 1; 1 0; 1 1];  % combinations
p_rk = prod(prob_gamma(Z', p_rk), 1)';
p_rk = p_rk ./ sum(p_rk);  % normalized
assert(isequal(round(p_rk, 6), ...
    [0.960596; 0.019504; 0.019504; 0.000396]))
assert(isequal(round(MKF_SF2.p_rk, 6), ...
    [0.960596; 0.019504; 0.019504; 0.000396]))

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5; -0.2; -0.4];
MKF_testx0 = MKFObserverSF_RODD(sys_model,io,P0,epsilon,sigma_wp, ...
    Q0,R,f,m,d,label,x0);
assert(isequal(MKF_testx0.xkp1_est, x0))

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
RMSE_results = containers.Map();
for i = 1:n_obs

    obs = observers{i};
    [obs, sim_result] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_result.X_est(1:5,:) - X(1:5, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_result.Y_est(1:5,:) - Y(1:5, :)) < 1e-10, [1 2]))

    % Check observer static errors are small
    % after input disturbance
    % TODO: Should these be closer?
    if all(sigma_MP == 0)
        assert(all(abs(sim_result.Y_est(end, :) - Y(end, :)) < 1e-3, [1 2]));
        assert(all(abs(sim_result.X_est(end, 3:4) - du0) < 1e-3, [1 2]));
    end

    % Compute mean-squared error
    Y_est = sim_result.Y_est;
    RMSE_results(obs.label) = mean((Y_est - Y).^2);
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

% Results on 2022-12-07 after fixing MKFObserverSF_DI
RMSE_test_values = containers.Map(...
 {'KF3',               'MKF_SF1',           'MKF_SF2', ...
  'MKF3',              'MKF4',              'SKF'}, ...
 {[0.000296 0.000433], [0.000852, 0.000739], [0.000291, 0.000328], ...
  [0.000971 0.000916], [0.000017 0.000022], [0.000017 0.000022]} ...
);

% for label = RMSE_results.keys
%     fprintf("%s: %f, %f (%f, %f)\n", label{1}, RMSE_results(label{1}), RMSE_test_values(label{1}))
% end
for label = RMSE_results.keys
    assert(isequal(round(RMSE_results(label{1}), 6), RMSE_test_values(label{1})))
end
