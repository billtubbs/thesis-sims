% Test GPB1_update.m

clear all
rng(1)


%% Test calculation - SISO system

% Define system

% Sample period
Ts = 0.5;

% Discrete time state space models
% Model #1
A1 = 0.7;
B1 = 1;
C1 = 0.3;
D1 = 0;
Gpss1 = ss(A1,B1,C1,D1,Ts);

% Model #2
A2 = 0.9;
B2 = 1;
C2 = -0.3;  % -ve gain!
D2 = 0;
Gpss2 = ss(A2,B2,C2,D2,Ts);

% Dimensions
n = size(A1, 1);
nu = size(B1, 2);
ny = size(C1, 1);

% Check dimensions
assert(isequal(size(A1), size(A2)))
assert(isequal(size(B1), size(B2)))
assert(isequal(size(C1), size(C2)))
assert(isequal(size(D1), size(D2)))

% Define system models
A = {A1, A2};
B = {B1, B2};
C = {C1, C2};
D = {D1, D2};
nj = numel(A);

% Input disturbance variance
%sigma_w = 0.1;
sigma_w = 0;

% Process noise std. dev.
sigma_W = [0; 0];

% Measurement noise std. dev.
sigma_M = 0.1;

% Transition probabilities
T = [0.95 0.05; 0.01 0.99];
assert(all(sum(T, 2) == 1))

% Observer parameters (same for all observers)
P0 = 10000;
x0 = 0.5;
Q1 = 0.01;
R1 = 0.1^2;
Q2 = 0.01;
R2 = 0.1^2;

% Define model structs
m1.A = A1;
m1.B = B1;
m1.C = C1;
m1.Q = Q1;
m1.R = R1;
m2.A = A2;
m2.B = B2;
m2.C = C2;
m2.Q = Q2;
m2.R = R2;
models = {m1, m2};

% Test GPB1 update and prediction steps

% Inputs
xk_est = x0;
uk = 0;
yk = 0.1830;
p_seq_g_Yk = [0.4; 0.6];
Pk = P0;
n_filt = nj;

% Do prediction step first to calculate priors
% Prepare empty arrays
Xkp1f_est = nan(n, 1, n_filt);
Pkp1f = nan(n, n, n_filt);
Ykp1f_est = nan(ny, 1, n_filt);
for j = 1:n_filt
    [Xkp1f_est(:,:,j), Pkp1f(:,:,j)] = ...
        kalman_predict_f( ...
            models{j}.A, models{j}.B, models{j}.Q, ...
            xk_est, Pk, uk ...
        );
    Ykp1f_est(:,:,j) = models{j}.C * Xkp1f_est(:,:,j);
end

% Update step
[xk_est1,yk_est,Pk1,p_seq_g_Yk1] = GPB1_update(models,T,Xkp1f_est, ...
    Pkp1f,yk,p_seq_g_Yk);

% Compare to code from Zhao et al.
% Note: original code does not include known inputs u(k)
% Also, D matrix is used for different purpose.
Model = struct;
Model.A = {A1, A2};
Model.Bu = {B1, B2};
Model.Bw = {1, 1};
Model.C = {C1, C2};
Model.D = repmat({eye(ny)},1,nj);
Model.Q = {Q1, Q2};
Model.R = {R1, R2};
Model.TP = T;
[x_test,P_test,mu_test] = GPB1_estimation(xk_est,Pk,yk,uk,Model,p_seq_g_Yk);

% Compare GPB1 estimates
assert(isequal(x_test, xk_est1))
assert(isequal(P_test, Pk1))
assert(isequal(mu_test, p_seq_g_Yk1))


% TODO: Finish GPB2 tests

% Test GPB2 update and prediction steps

% Inputs
xk_est = x0;
uk = 0;
yk = 0.1830;

% Merged estimates from last time instant
xki_est = repmat(x0, 1, 1, nj);
xki_est(:,:,2) = xki_est(:,:,2) + 0.1;  % make a slight difference in xk2_est
Pki = repmat(P0, 1, 1, nj);
Pki(:,:,2) = Pki(:,:,2) + 1000;  % slight difference in xk2_est
p_seqi_g_Yk = [0.4; 0.6];

% Mode transitions
gamma_km1 = [0 1 0 1];
gamma_k = [0 0 1 1];

% Do prediction step first to calculate prior estimates (k|k-1)
n_filt = nj*nj;
for f = 1:n_filt
    i = gamma_km1(f) + 1;
    j = gamma_k(f) + 1;
    m = models{j};  % TODO: Change to gamma_km1 later
    [Xkp1f_est(:,:,j), Pkp1f(:,:,j)] = ...
        kalman_predict_f( ...
            models{j}.A, models{j}.B, models{j}.Q, ...
            xki_est(:,:,i), Pki(:,:,i), uk ...
        );
    Ykp1f_est(:,:,j) = models{j}.C * Xkp1f_est(:,:,j);
end

% Update step
% TODO: This is not correct - check GPB2_update after adding uk.
% [xk_est_out,yk_est_out,Pk_out,xk_est2,yk_est2,Pk2,p_seq_g_Yk2] = ...
%       GPB2_update(models,T,Xkp1f_est,Pkp1f,yk,p_seq_g_Yk);
% 
% % Compare to code from Zhao et al.
% % Note: original code does not include known inputs u(k)
% x = xki_est;
% P = Pki;
% [x_test,P_test,mu_test,Out_x_test,Out_P_test] = ...
%     GPB2_estimation(x,P,yk,uk,Model,p_seq_g_Yk);
% 
% assert(isequal(round(x_test, 4), [0.61 -0.61]))
% assert(isequal(round(P_test, 4), cat(3, 0.1111, 0.1111)))
% assert(isequal(round(mu_test, 4), [0.4582 0.5418]))
% assert(isequal(round(Out_x_test, 6), -0.050945))
% assert(isequal(round(Out_P_test, 6), 0.480601))
% 
% % Compare GPB2 estimates
% assert(isequal(x_test, xk_est2))
% assert(isequal(P_test, Pk2))
% assert(isequal(u_test, p_seq_g_Yk2))
% assert(isequal(Out_x_test, xk_est_out))
% assert(isequal(Out_P_test, Pk1))
% 




%% Test calculation - 2x2 system

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
C = [ 0.1110 0         0  0;
             0  0.1110 0  0];
D = zeros(2, 4);

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Designate measured input and output signals
u_meas = [true; true; false; false];
y_meas = [true; true];

% Disturbance inputs to system
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.1; 0.1];
sigma_wp = [0.01 1; 0.01 1];

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);
model = struct;
model.A = A;
model.B = Bu;
model.C = C;
assert(all(D == 0, [1 2]))
nj = 3;
models = repmat({model}, 1, nj);
models{1}.Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]);
models{2}.Q = diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]);
models{3}.Q = diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2]);
models{1}.R = diag(sigma_M.^2); 
models{2}.R = diag(sigma_M.^2);
models{3}.R = diag(sigma_M.^2);

% Observer initialization parameters
x0 = [0.5; -0.1; 0.05; -0.05];
P0 = 1000*eye(n);

p_gamma = [1-epsilon epsilon]';
Z = [0 0; 1 0; 0 1];  % combinations
p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
p_gamma = p_gamma ./ sum(p_gamma);  % normalized
T = repmat(p_gamma', 3, 1);

% Compare to code from Zhao et al.
Model = struct;
Model.A = cellfun(@(m) m.A, models, 'UniformOutput', false);
Model.Bu = cellfun(@(m) m.B, models, 'UniformOutput', false);
Model.Bw = repmat({eye(n)},1,nj);  % Is this correct?
Model.C = cellfun(@(m) m.C, models, 'UniformOutput', false);
Model.D = repmat({eye(ny)},1,nj);  % D is used for different purpose
Model.Q = cellfun(@(m) m.Q, models, 'UniformOutput', false);
Model.R = cellfun(@(m) m.R, models, 'UniformOutput', false);
Model.TP = T;
xk_est = x0;
yk = sigma_M .* randn(ny,1);
p_seq_g_Yk = [0.6; 0.2; 0.2];
Pk = P0;
uk = zeros(nu,1);
[x_test,P_test,mu_test] = GPB1_estimation(xk_est,Pk,yk,uk,Model,p_seq_g_Yk);

% Prepare cell array of structs to store data
nh = nj;
Xkp1f_est = nan(n, 1, nh);
Pkp1f = nan(n, n, nh);
Ykp1f_est = nan(ny, 1, nh);

% Do prediction step first
for j = 1:nh
    [Xkp1f_est(:,:,j), Pkp1f(:,:,j)] = ...
        kalman_predict_f( ...
            models{j}.A, models{j}.B, models{j}.Q, ...
            xk_est, Pk, uk ...
        );
    Ykp1f_est(:,:,j) = models{j}.C * Xkp1f_est(:,:,j);
end

% Update step
[xk_est,yk_est,Pk,p_seq_g_Yk] = GPB1_update(models,T,Xkp1f_est, ...
    Pkp1f,yk,p_seq_g_Yk);

% Compare
assert(isequal(x_test, xk_est))
assert(max(abs(P_test - Pk), [], "all") < 1e-16)
assert(isequal(mu_test, p_seq_g_Yk))