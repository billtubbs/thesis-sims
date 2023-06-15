% test get_obs_params.m

clear all


%% Test all fieldnames returned

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

params = get_obs_params(KFPSS1);
assert(isequal(fieldnames(params), {'model', 'K'}'))
assert(isequal(params.model, KFPSS1.model))
assert(isequal(params.K, KFPSS1.K))

params = get_obs_params(KFFSS1);
assert(isequal(fieldnames(params), {'model', 'Kf'}'))
assert(isequal(params.model, KFFSS1.model))
assert(isequal(params.Kf, KFFSS1.Kf))

params = get_obs_params(KF1);
assert(isequal(fieldnames(params), {'model', 'P0'}'))
assert(isequal(params.model, KF1.model))
assert(isequal(params.P0, KF1.P0))

params = get_obs_params(LB1);
assert(isequal(fieldnames(params), {'model', 'poles', 'K'}'))
assert(isequal(params.model, LB1.model))
assert(isequal(params.poles, LB1.poles))
assert(isequal(params.K, LB1.K))

% Define a switching Kalman filter
P0 = 1000*eye(n);
Q1 = diag([Q0(1,1) sigma_wp{1}(1)^2]);
Q2 = diag([Q0(1,1) sigma_wp{1}(2)^2]);
R = sigma_M^2;
models = {model, model};
models{1}.Bu = Bu;
models{2}.Bu = Bu;
models{1}.Q = Q1;
models{2}.Q = Q2;
models{1}.R = R;
models{2}.R = R;
SKF1 = SKFObserver(models,P0,"SKF1");
params = get_obs_params(SKF1);
assert(isequal(fieldnames(params), {'models', 'P0', 'nj'}'))
assert(isequal(params.models, SKF1.models))
assert(isequal(params.P0, SKF1.P0))
assert(isequal(params.nj, SKF1.nj))

% Simulation parameters
nT = 100;
t = Ts*(0:nT)';

% Define a switching Kalman filter with sequence
seq = {zeros(1, nT+1)};
seq{1}(t == 10) = 1;
SKF2 = SKFObserverS(models,P0,seq{1},"SKF2");
params = get_obs_params(SKF2);
assert(isequal(fieldnames(params), {'models', 'P0', 'nj'}'))
assert(isequal(params.models, SKF2.models))
assert(isequal(params.P0, SKF2.P0))
assert(isequal(params.nj, SKF2.nj))

params = get_obs_params(MKF_SF95);
assert(isequal(fieldnames(params), {'model', 'P0', 'Q0', 'R', ...
    'epsilon', 'sigma_wp', 'nf', 'f', 'm', 'd', 'nh', 'nm', ...
    'nh_max', 'beta'}'))
assert(isequal(params.P0, MKF_SF95.P0))
assert(isequal(params.Q0, MKF_SF95.Q0))
assert(isequal(params.R, MKF_SF95.R))
assert(isequal(params.epsilon, MKF_SF95.epsilon))
assert(isequal(params.sigma_wp, MKF_SF95.sigma_wp))
assert(isequal(params.nf, MKF_SF95.nf))
assert(isequal(params.f, MKF_SF95.f))
assert(isequal(params.m, MKF_SF95.m))
assert(isequal(params.d, MKF_SF95.d))
assert(isequal(params.nh, MKF_SF95.nh))
assert(isequal(params.nm, MKF_SF95.nm))
assert(isequal(params.nh_max, MKF_SF95.nh_max))
assert(isequal(params.beta, MKF_SF95.beta))

params = get_obs_params(MKF_SF1);
assert(isequal(fieldnames(params), {'model', 'P0', 'Q0', 'R', ...
    'epsilon', 'sigma_wp', 'nf', 'f', 'm', 'd', 'nh', 'nm', ...
    'nh_max', 'beta'}'))
assert(isequal(params.model, MKF_SF1.sys_model))
assert(isequal(params.P0, MKF_SF1.P0))
assert(isequal(params.Q0, MKF_SF1.Q0))
assert(isequal(params.R, MKF_SF1.R))
assert(isequal(params.epsilon, MKF_SF1.epsilon))
assert(isequal(params.sigma_wp, MKF_SF1.sigma_wp))
assert(isequal(params.nf, MKF_SF1.nf))
assert(isequal(params.f, MKF_SF1.f))
assert(isequal(params.m, MKF_SF1.m))
assert(isequal(params.d, MKF_SF1.d))
assert(isequal(params.nh, MKF_SF1.nh))
assert(isequal(params.nm, MKF_SF1.nm))
assert(isequal(params.nh_max, MKF_SF1.nh_max))
assert(isequal(params.beta, MKF_SF1.beta))

params = get_obs_params(MKF_SP1);
assert(isequal(fieldnames(params), {'model', 'P0', 'Q0', 'R', 'epsilon', ...
    'sigma_wp', 'nh', 'n_min'}'))
assert(isequal(params.model, MKF_SP1.sys_model))
assert(isequal(params.P0, MKF_SP1.P0))
assert(isequal(params.Q0, MKF_SP1.Q0))
assert(isequal(params.R, MKF_SP1.R))
assert(isequal(params.epsilon, MKF_SP1.epsilon))
assert(isequal(params.sigma_wp, MKF_SP1.sigma_wp))
assert(isequal(params.nh, MKF_SP1.nh))
assert(isequal(params.n_min, MKF_SP1.n_min))

% Multiple model observer with sequence pruning
nh = 10;  % number of filters
n_min = 7;  % minimum life of cloned filters
epsilon = 0.05;
T = [1-epsilon epsilon; epsilon 1-epsilon];
MKF_SP = MKFObserverSP(models,P0,T,nh,n_min,"MKF_SP");

params = get_obs_params(MKF_SP);
assert(isequal(fieldnames(params), {'models', 'P0', 'T', 'nh', 'n_min'}'))
assert(isequal(params.models, MKF_SP.models))
assert(isequal(params.T, MKF_SP.T))
assert(isequal(params.nh, MKF_SP.nh))
assert(isequal(params.n_min, MKF_SP.n_min))

% Multiple model observer with sequence pruning
% and detection interval
d = 1;
MKF_SP_DI = MKFObserverSP_DI(models,P0,T,d,nh,n_min,"MKF_SP_DI");

params = get_obs_params(MKF_SP_DI);
assert(isequal(fieldnames(params), {'models', 'P0', 'T', 'd', 'nh', ...
    'n_min'}'))
assert(isequal(params.models, MKF_SP_DI.models))
assert(isequal(params.T, MKF_SP_DI.T))
assert(isequal(params.d, MKF_SP_DI.d))
assert(isequal(params.nh, MKF_SP_DI.nh))
assert(isequal(params.n_min, MKF_SP_DI.n_min))


%% Test other MKF observers

% Load switching system definition
sys_js2_siso

% Transition probabilities
T = [0.95 0.05; 0.01 0.99];
assert(all(sum(T, 2) == 1))

% Observer parameters
P0 = 10000;
models{1}.Q = 0.01;
models{1}.R = 0.1^2;
models{2}.Q = 0.01;
models{2}.R = 0.1^2;

% AMM
MKF_AMM1 = MKFObserverAMM(models,P0,"MKF_AMM1",x0);
params = get_obs_params(MKF_AMM1);
assert(isequal(fieldnames(params), {'models', 'P0', 'nh'}'))
assert(isequal(params.models, MKF_AMM1.models))
assert(isequal(params.nh, MKF_AMM1.nh))

% GPB1
MKF_GPB1 = MKFObserverGPB1(models,P0,T,"MKF_GPB1",x0);
params = get_obs_params(MKF_GPB1);
assert(isequal(fieldnames(params), {'models', 'P0', 'T', 'nh'}'))
assert(isequal(params.models, MKF_GPB1.models))
assert(isequal(params.T, MKF_GPB1.T))
assert(isequal(params.nh, MKF_GPB1.nh))

% GPB2
MKF_GPB2 = MKFObserverGPB2(models,P0,T,"MKF_GPB2",x0);
params = get_obs_params(MKF_GPB2);
assert(isequal(fieldnames(params), {'models', 'P0', 'T', 'nh'}'))
assert(isequal(params.models, MKF_GPB2.models))
assert(isequal(params.T, MKF_GPB2.T))
assert(isequal(params.nh, MKF_GPB2.nh))

