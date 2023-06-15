% Tests the following classes:
% 
%  - MKFObserverSP
%  - MKFObserverSP_DI
%  - MKFObserverSP_RODD
%

clear all

% If plotting:
%plot_dir = 'plots';
%addpath("../plot-utils")

seed = 0;
rng(seed)


%% Test initialization with SISO system

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Also define two MKF_SP observers with detection intervals
T = repmat([1-epsilon epsilon], 2, 1);
d = 1;  % detection interval — this should produce identical results to MKF_SP
nh = 10;
n_min = 7;
MKF_SP_DI1 = MKFObserverSP_DI(obs_models,P0,T,d,nh,n_min,"MKF_SP_DI1");

d = 5;  % detection interval
nh = 10;
n_min = 7;
MKF_SP_DI2 = MKFObserverSP_DI(obs_models,P0,T,d,nh,n_min,"MKF_SP_DI2");

% Check observer attributes
assert(strcmp(MKF_SP_DI1.type, "MKF_SP_DI"))
assert(isequal(MKF_SP_DI1.T, T))
assert(MKF_SP_DI1.d == 1)
assert(MKF_SP_DI1.nh == 10)
assert(MKF_SP_DI1.n_min == 7)
assert(isequal(MKF_SP_DI1.n_hold, 7))
assert(isequal(MKF_SP_DI1.n_main, 3))
assert(isequaln(MKF_SP_DI1.f_hold, 4:10))
assert(isequaln(MKF_SP_DI1.f_main, 1:3))
assert(isequaln(MKF_SP_DI1.id, 0))
assert(isequaln(MKF_SP_DI1.id_next, 1))
assert(MKF_SP_DI1.n == 2)
assert(MKF_SP_DI1.nu == 1)
assert(MKF_SP_DI1.ny == 1)
assert(MKF_SP_DI1.nj == 2)
assert(isequal(MKF_SP_DI1.models, obs_models))
assert(isequal(MKF_SP_DI1.Ts, Ts))
assert(isequal(MKF_SP_DI1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SP_DI1.Pkp1, 1000*eye(2)))
assert(isequal(MKF_SP_DI1.r0, ones(MKF_SP_DI1.nh, 1)))
assert(isequal(MKF_SP_DI1.p_seq_g_Yk_init, ones(MKF_SP_DI1.nh, 1)./nh))
assert(isequaln(MKF_SP_DI1.p_rk_g_rkm1, nan(MKF_SP_DI1.nh, 1)))

% Don't bother checking every attribute of this one
assert(MKF_SP_DI2.d == 5)
% TODO: Include above in simulations

% Observer model without disturbance noise input
Bu = B(:, u_known);
Du = D(:, u_known);
nu = sum(u_known);
nw = sum(~u_known);

% Set noise variances for observer design
sigma_M = 0.1;
sigma_W = [0; 0];

% Check observer attributes
assert(strcmp(MKF_SP1.type, "MKF_SP_RODD"))
assert(MKF_SP1.epsilon == 0.01)
assert(isequal(MKF_SP1.sigma_wp, sigma_wp))
assert(isequal(MKF_SP1.Q0, Q0))
assert(isequal(MKF_SP1.R, R))
assert(MKF_SP1.nh == 10)
assert(MKF_SP1.n_min == 7)
assert(isequal(MKF_SP1.n_hold, 7))
assert(isequal(MKF_SP1.n_main, 3))
assert(isequaln(MKF_SP1.f_hold, 4:10))
assert(isequaln(MKF_SP1.f_main, 1:3))
% assert(isequaln(MKF_SP1.i, 0))
assert(MKF_SP1.n == 2)
assert(MKF_SP1.nu == 1)
assert(MKF_SP1.ny == 1)
assert(MKF_SP1.nj == 2)
assert(isequal(MKF_SP1.sys_model, model))
assert(isequal(MKF_SP1.Ts, Ts))
assert(isequaln(MKF_SP1.io.u_known, u_known))
assert(isequal(MKF_SP1.models{1}.Q, [0.01 0; 0 sigma_wp{1}(1)^2]))
assert(isequal(MKF_SP1.models{2}.Q, [0.01 0; 0 sigma_wp{1}(2)^2]))
assert(isequal(MKF_SP1.models{1}.R, R))
assert(isequal(MKF_SP1.models{2}.R, R))
% assert(isequal(size(MKF_SP1.seq), [MKF_SP1.nh 1]))
% assert(isequal(size(cell2mat(MKF_SP1.seq)), [MKF_SP1.nh MKF_SP1.nf]))
% assert(MKF_SP1.nf == size(MKF_SP1.seq{1}, 2))
assert(isequal(MKF_SP1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SP1.Pkp1, 1000*eye(2)))
assert(isequal(MKF_SP1.r0, ones(MKF_SP1.nh, 1)))
assert(isequal(MKF_SP1.p_seq_g_Yk_init, [1; zeros(MKF_SP1.nh-1, 1)]))
assert(isequaln(MKF_SP1.p_rk_g_rkm1, nan(MKF_SP1.nh, 1)))

% Check initialization of filters
assert(isequal(MKF_SP1.filters.Pkp1(:,:,1), MKF_SP1.P0))
assert(isequal(MKF_SP1.filters.Xkp1_est(:,:,1), MKF_SP1.x0))
for i = 2:MKF_SP1.nh
    assert(isequal(MKF_SP1.filters.Pkp1(:,:,i), 1e10*eye(2)))
    assert(isequal(MKF_SP1.filters.Xkp1_est(:,:,i), MKF_SP1.x0))
end

% Check observer attributes
assert(strcmp(MKF_SP2.type, "MKF_SP_RODD"))
assert(MKF_SP2.epsilon == 0.01)
assert(isequal(MKF_SP2.sigma_wp, sigma_wp))
assert(isequal(MKF_SP2.Q0, Q0))
assert(isequal(MKF_SP2.R, R))
assert(MKF_SP2.nh == 25)
assert(MKF_SP2.n_min == 21)
assert(isequal(MKF_SP2.n_hold, 21))
assert(isequal(MKF_SP2.n_main, 4))
assert(isequaln(MKF_SP2.f_hold, 5:25))
assert(isequaln(MKF_SP2.f_main, 1:4))
% assert(isequaln(MKF_SP2.i, 0))
assert(MKF_SP2.n == 2)
assert(MKF_SP2.nu == 1)
assert(MKF_SP2.ny == 1)
assert(MKF_SP2.nj == 2)
assert(isequal(MKF_SP2.sys_model, model))
assert(isequal(MKF_SP2.Ts, Ts))
assert(isequaln(MKF_SP2.io.u_known, u_known))
assert(isequal(MKF_SP2.models{1}.Q, [0.01 0; 0 sigma_wp{1}(1)^2]))
assert(isequal(MKF_SP2.models{2}.Q, [0.01 0; 0 sigma_wp{1}(2)^2]))
assert(isequal(MKF_SP2.models{1}.R, R))
assert(isequal(MKF_SP2.models{2}.R, R))
% assert(isequal(size(MKF_SP2.seq), [MKF_SP2.nh 1]))
% assert(isequal(size(cell2mat(MKF_SP2.seq)), [MKF_SP2.nh MKF_SP2.f]))
% assert(MKF_SP2.nf == size(MKF_SP2.seq{1}, 2))
assert(isequal(MKF_SP2.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SP2.Pkp1, 1000*eye(2)))
assert(isequal(MKF_SP2.r0, ones(MKF_SP2.nh, 1)))
assert(isequal(MKF_SP2.p_seq_g_Yk_init, [1; zeros(MKF_SP2.nh-1, 1)]))
assert(isequaln(MKF_SP2.p_rk_g_rkm1, nan(MKF_SP2.nh, 1)))

% Check initialization of filters
assert(isequal(MKF_SP1.filters.Pkp1(:,:,1), MKF_SP1.P0))
assert(isequal(MKF_SP1.filters.Xkp1_est(:,:,1), MKF_SP1.x0))
for i = 2:MKF_SP1.nh
    assert(isequal(MKF_SP1.filters.Pkp1(:,:,i), 1e10*eye(2)))
    assert(isequal(MKF_SP1.filters.Xkp1_est(:,:,i), MKF_SP1.x0))
end

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5];
io.u_known = u_known;
io.y_meas = true(ny, 1);
d = 1;
MKF_SP_testx0 = MKFObserverSP_RODD(model,io,P0,epsilon, ...
                sigma_wp,Q0,R,nh,n_min,label,d,x0);
assert(isequal(MKF_SP_testx0.xkp1_est, x0))
assert(isequal( ...
    MKF_SP_testx0.filters.Xkp1_est, repmat(x0, 1, 1, MKF_SP_testx0.nh) ...
))
assert(isequal(MKF_SP_testx0.r0, ones(MKF_SP_testx0.nh, 1)))
assert(isequal(MKF_SP_testx0.p_seq_g_Yk_init, [1; zeros(MKF_SP_testx0.nh-1, 1)]))


%% Test convergence to steady-state

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Check steady-state at x0 = [0; 0]
obs = MKF_SP1;
assert(isequal(obs.xkp1_est, [0; 0]))
nT = 10;
U_m = zeros(nT+1, sum(u_known));
Y_m = zeros(nT+1, ny);
for i = 1:(nT+1)
    uk = U_m(i,:)';
    yk = Y_m(i,:)';
    obs.update(yk, uk);
    assert(isequal(obs.xk_est, [0; 0]))
    assert(isequal(obs.yk_est, 0))
end

% Check steady-state at x0 = [1; 0]
x0 = [1; 0];
io.u_known = u_known;
io.y_meas = true(ny, 1);
d = 1;
obs = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nh,n_min, ...
    label,d,x0);
assert(isequal(obs.xkp1_est, x0))
nT = 10;
U_m = 0.3*ones(nT+1, 1);
U = [U_m zeros(nT+1,1)];
t = Ts*(0:nT)';
[Y,~,X] = lsim(Gpss,U,t,x0);
assert(all(Y == Y(1,1)))
for i = 1:(nT+1)
    uk = U_m(i,:)';
    yk = Y(i,:)';
    obs.update(yk, uk);
    assert(all(abs(obs.xk_est - x0) < 1e-6))
    assert(abs(obs.yk_est - 0.3) < 1e-6)
end


%% Test sequence updates on SISO linear system

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

x0 = [0; 0];
io.u_known = u_known;
io.y_meas = true(ny, 1);
nh = 5;
n_min = 2;
d = 1;
obs = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nh,n_min, ...
    label,d,x0);
assert(isequal(obs.xkp1_est, x0))

% % Generate test simulation data
% nT = 10;
% U_m = zeros(nT+1, sum(u_known));
% % Add a random shock
% Wp = zeros(nT+1, sum(~u_known));
% Wp(5, :) = 1;
% % Compute outputs (n0 measurement noise)
% [Y_m, t] = lsim(Gpss, [U_m Wp], t);
% [t U_m Wp Y_m]

% Test simuation data
sim_data = [ ...
         0         0         0         0;
    0.5000         0         0         0;
    1.0000         0         0         0;
    1.5000         0         0         0;
    2.0000         0    1.0000         0;
    2.5000         0         0         0;
    3.0000         0         0    0.3000;
    3.5000         0         0    0.5100;
    4.0000         0         0    0.6570;
    4.5000         0         0    0.7599;
    5.0000         0         0    0.8319];
nT = size(sim_data, 1);
t = sim_data(:, 1);
U_m = sim_data(:, 2);
Wp = sim_data(:, 3);
Y_m = sim_data(:, 4);

% Set marker values on each sequence - for testing only
% these values at the end of the sequences are not used
% by the observer.
% for i = 1:nh
%     obs.seq{i}(8) = i;
% end
seq0 = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    0 0 0 0 0 0 0 3
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 5
];
%%disp(obs.i)  % use for debugging
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 0))
% assert(isequaln(obs.i_next, 1))
% assert(isequaln(cell2mat(obs.seq), seq0))
assert(isequal(obs.n_hold, 2))
assert(isequal(obs.n_main, 3))
assert(isequaln(obs.f_hold, [4 5]))
assert(isequaln(obs.f_main, [1 2 3]))

% Check initialization of probabilities
assert(isequal(obs.rk, [1 1 1 1 1]'))
assert(isequaln(obs.p_rk_g_rkm1, nan(5, 1)))
assert(isequaln(obs.p_yk_g_seq_Ykm1, nan(5, 1)))
assert(isequaln(obs.p_seq_g_Ykm1, nan(5, 1)))
assert(isequal(obs.p_seq_g_Yk, [1 zeros(1, 4)]'))

% Check initialization of filters
assert(isequal(obs.filters.Pkp1(:,:,1), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,2), 1e10*eye(2)))
assert(isequal(obs.filters.Pkp1(:,:,3), 1e10*eye(2)))
assert(isequal(obs.filters.Pkp1(:,:,4), 1e10*eye(2)))
assert(isequal(obs.filters.Pkp1(:,:,5), 1e10*eye(2)))

% Check initial estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, nan(2, 1)))
assert(isequaln(obs.yk_est, nan))

% Update at k = 0
i = 1;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    1 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 5
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 1))
% assert(isequaln(obs.i_next, 2))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [5 3]))
assert(isequaln(obs.f_main, [1 2 4]))

% Check probabilities
assert(isequal(obs.rk, [1 1 2 1 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.042 0 0.042 0 0]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.99 0 0.01 0 0]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.99 0 0.01 0 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 1
i = 2;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    1 0 0 0 0 0 0 1
    0 1 0 0 0 0 0 1
    0 0 0 0 0 0 0 5
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 2))
% assert(isequaln(obs.i_next, 3))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [3 4]))
assert(isequaln(obs.f_main, [1 2 5]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 2 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.01 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0420 0 0.0420 0.0420 0]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9801 0 0.0099 0.0099 0]'))  % NOTE: doesn't quite sum to 1
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9802 0 0.0099 0.0099 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 2
i = 3;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    1 0 0 0 0 0 0 1
    0 1 0 0 0 0 0 1
    0 0 1 0 0 0 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 3))
% assert(isequaln(obs.i_next, 4))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [4 5]))
assert(isequaln(obs.f_main, [1 2 3]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 2]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [1.8682 1.868 1.0834 1.8682 1.8682]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9704 0 0.0098 0.0098 0.0098]'))  % NOTE: doesn't quite sum to 1
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9746 0 0.0057 0.0098 0.0098]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 3
i = 4;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1
    1 0 0 0 0 0 0 1
    0 1 0 0 0 0 0 1
    0 0 1 0 0 0 0 1  % all sequences are now splits from #1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 4))
% assert(isequaln(obs.i_next, 5))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [5 2]))
assert(isequaln(obs.f_main, [1 3 4]))

% Check probabilities
assert(isequal(obs.rk, [1 2 1 1 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.01 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [2.4676 2.4676 2.0186 1.1707 2.4676]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9649 0.0097 0.0057 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9711 0.0098 0.0047 0.0047 0.0098]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 4
i = 5;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 0 0 0 1 0 0 1  % this seq. has now been replaced
    0 1 0 0 0 0 0 1
    0 0 1 0 0 0 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 5))
% assert(isequaln(obs.i_next, 6))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [2 3]))
assert(isequaln(obs.f_main, [1 4 5]))

% Check probabilities
assert(isequal(obs.rk, [1 1 2 1 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [2.8078 2.8078 2.8078 2.0187 1.2019]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9614 0.0097 0.0097 0.0046 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9728 0.0098 0.0098 0.0034 0.0042]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 5
i = 6;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 0 0 0 1 0 0 1
    0 0 0 0 0 1 0 1
    0 0 1 0 0 0 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 6))
% assert(isequaln(obs.i_next, 7))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [3 4]))
assert(isequaln(obs.f_main, [1 5 2]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 2 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.01 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [3.0218 1.2172 3.0218 3.0218 2.0223]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9631 0.0097 0.0097 0.0097 0.0042]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9736 0.004 0.0098 0.0098 0.0028]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 6  *** First non-zero measurement ***
i = 7;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 0 0 0 1 0 0 1
    0 0 0 0 0 1 0 1
    0 0 0 0 0 0 1 1
];
%disp(obs.i)
%display_obs_data(obs)  
%disp(debug_array(obs))
% assert(isequaln(obs.i, 7))
% assert(isequaln(obs.i_next, 8))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [4 5]))
assert(isequaln(obs.f_main, [1 2 3]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 2]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.1865 0.6345 0.8015 0.1865 0.1865]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9638 0.0039 0.0097 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9281 0.0129 0.0403 0.0094 0.0094]'))

% For comparison: probability densities if yk had been 0:
% assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [3.1646 1.2260 3.1646 3.1646 3.1646]'))
% assert(isequal(round(obs.p_seq_g_Yk, 4), [0.0050    0.0019    0.0050    0.4940    0.4940]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [0.4290; 0.1510]))
assert(isequaln(round(obs.xk_est, 4), [0.3971; 0.1510]))
assert(isequaln(round(obs.yk_est, 4), 0.1191))
assert(isequal(round(obs.yk_est - yk, 4), -0.1809))

% Update at k = 7
i = 8;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1
    0 0 0 0 1 0 0 0
    0 0 0 0 0 1 0 0
    0 0 0 0 0 0 1 0
];
%disp(obs.i)
%display_obs_data(obs)  
%disp(debug_array(obs))
% assert(isequaln(obs.i, 8))
% assert(isequaln(obs.i_next, 1))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [5 2]))
assert(isequaln(obs.f_main, [1 3 4]))

% Check probabilities
assert(isequal(obs.rk, [1 2 1 1 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.01 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0166 0.0166 1.9388 0.5807 0.0166]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9188 0.0093 0.0399 0.0093 0.0093]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.1553 0.0016 0.7867 0.0548 0.0016]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [1.9184; 0.8599]))
assert(isequaln(round(obs.xk_est, 4), [1.5122; 0.8599]))
assert(isequaln(round(obs.yk_est, 4), 0.4537))
assert(isequal(round(obs.yk_est - yk, 4), -0.0563))

% Update at k = 8
i = 9;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1
    0 0 0 0 1 0 0 0
    1 0 0 0 1 0 0 0  % New additions have looped back to 1st position
    0 0 0 0 0 0 1 0
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 1))
% assert(isequaln(obs.i_next, 2))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [2 4]))
assert(isequaln(obs.f_main, [3 1 5]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 2 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.01 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0084 0.0084 2.5317 2.5317 0.5437]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.1537 0.0016 0.7789 0.0079 0.0016]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.0006 0 0.9889 0.01 0.0004]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [2.4869; 0.9775]))
assert(isequaln(round(obs.xk_est, 4), [2.1563; 0.9775]))
assert(isequaln(round(obs.yk_est, 4), 0.6469))
assert(isequal(round(obs.yk_est - yk, 4), -0.0101))


%% Test sequence updates 2

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

x0 = [0; 0];
nh = 5;
n_min = 1;  % NOTE: this produces identical results to previous
            % MKFObserverSP_RODD and mkf_observer_AFMM with n_min = 2;
            % This is due to the interpretation about whether a
            % hypothesis leaving holding group goes into main group
            % first (as in this version) or can be immediately 
            % eliminated before going to main group (as previously).
io.u_known = u_known;
io.y_meas = true(ny, 1);
d = 1;
obs = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nh, ...
    n_min,label,d,x0);
assert(isequal(obs.xkp1_est, x0))

% % Generate test simulation data
% nT = 10;
% U_m = zeros(nT+1, sum(u_known));
% % Add a random shock
% Wp = zeros(nT+1, sum(~u_known));
% Wp(5, :) = 1;
% % Compute outputs (n0 measurement noise)
% [Y_m, t] = lsim(Gpss, [U_m Wp], t);
% [t U_m Wp Y_m]

% Test simuation data
sim_data = [ ...
         0         0         0         0;
    0.5000         0         0         0;
    1.0000         0         0         0;
    1.5000         0         0         0;
    2.0000         0    1.0000         0;
    2.5000         0         0         0;
    3.0000         0         0    0.3000;
    3.5000         0         0    0.5100;
    4.0000         0         0    0.6570;
    4.5000         0         0    0.7599;
    5.0000         0         0    0.8319];
nT = size(sim_data, 1);
t = sim_data(:, 1);
U_m = sim_data(:, 2);
Wp = sim_data(:, 3);
Y_m = sim_data(:, 4);

% Set marker values on each sequence - for testing only
% these values at the end of the sequences are not used
% by the observer.
% for i = 1:nh
%     obs.seq{i}(8) = i;
% end
seq0 = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    0 0 0 0 0 0 0 3
    0 0 0 0 0 0 0 4
    0 0 0 0 0 0 0 5
];
%%disp(obs.i)  % use for debugging
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 0))
% assert(isequaln(obs.i_next, 1))
% assert(isequaln(cell2mat(obs.seq), seq0))
assert(isequal(obs.n_hold, 1))
assert(isequal(obs.n_main, 4))
assert(isequaln(obs.f_hold, 5))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 1]'))
assert(isequaln(obs.p_rk_g_rkm1, nan(5, 1)))
assert(isequaln(obs.p_yk_g_seq_Ykm1, nan(5, 1)))
assert(isequaln(obs.p_seq_g_Ykm1, nan(5, 1)))
assert(isequal(obs.p_seq_g_Yk, [1 zeros(1, 4)]'))

% Check initialization of filters
assert(isequal(obs.filters.Pkp1(:,:,1), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,2), 1e10*eye(2)))
assert(isequal(obs.filters.Pkp1(:,:,3), 1e10*eye(2)))
assert(isequal(obs.filters.Pkp1(:,:,4), 1e10*eye(2)))
assert(isequal(obs.filters.Pkp1(:,:,5), 1e10*eye(2)))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, nan(2, 1)))
assert(isequaln(obs.yk_est, nan))

% Update at k = 0
i = 1;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    0 0 0 0 0 0 0 3
    1 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 5
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 1))
% assert(isequaln(obs.i_next, 2))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 4))
assert(isequaln(obs.f_main, [1 2 3 5]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 2 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.01 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.042 0 0 0.042 0]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.99 0 0 0.01 0]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.99 0 0 0.01 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequaln(obs.yk_est - yk, 0))

% Update at k = 1
i = 2;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    0 0 0 0 0 0 0 3
    1 0 0 0 0 0 0 1
    0 1 0 0 0 0 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 2))
% assert(isequaln(obs.i_next, 3))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 5))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 2]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.042 0 0 0.042 0.042]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9801 0 0 0.0099 0.0099]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9802 0 0 0.0099 0.0099]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequaln(obs.yk_est - yk, 0))

% Update at k = 2
i = 3;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 2
    0 0 1 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 1 0 0 0 0 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 3))
% assert(isequaln(obs.i_next, 4))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 3))
assert(isequaln(obs.f_main, [1 4 2 5]))

% Check probabilities
assert(isequal(obs.rk, [1 1 2 1 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [1.8682 1.868 1.8682 1.0834 1.8682]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9704 0 0.0098 0.0098 0.0098]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9746 0 0.0098 0.0057 0.0098]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequaln(obs.yk_est - yk, 0))

% Update at k = 3
i = 4;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 0 1 0 0 0 0 1
    1 0 0 0 0 0 0 1
    0 1 0 0 0 0 0 1  % all sequences are now splits from #1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 4))
% assert(isequaln(obs.i_next, 5))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 2))
assert(isequaln(obs.f_main, [1 5 4 3]))

% Check probabilities
assert(isequal(obs.rk, [1 2 1 1 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.01 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [2.4676 2.4676 2.4676 2.0186 1.1707]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9649 0.0097 0.0097 0.0057 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9711 0.0098 0.0098 0.0047 0.0047]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequaln(obs.yk_est - yk, 0))

% Update at k = 4
i = 5;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1  % this seq. has now been replaced
    0 0 1 0 0 0 0 1
    0 0 0 0 1 0 0 1
    0 1 0 0 0 0 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 5))
% assert(isequaln(obs.i_next, 6))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 4))
assert(isequaln(obs.f_main, [1 3 5 2]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 2 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.01 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [2.8078 2.8078 1.2019 2.8078 2.0187]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9614 0.0097 0.0097 0.0097 0.0046]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9728 0.0098 0.0042 0.0098 0.0034]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequaln(obs.yk_est - yk, 0))

% Update at k = 5
i = 6;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 0 1 0 0 0 0 1
    0 0 0 0 1 0 0 1
    0 0 0 0 0 1 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 6))
% assert(isequaln(obs.i_next, 7))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 5))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 2]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [3.0218 1.2172 2.0223 3.0218 3.0218]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9631 0.0097 0.0042 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9736 0.0040 0.0028 0.0098 0.0098]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0; 0]))
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequaln(obs.yk_est - yk, 0))

% Update at k = 6  *** First non-zero measurement ***
i = 7;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 1
    0 0 0 1 0 0 0 1
    0 0 0 0 0 0 1 1
    0 0 0 0 1 0 0 1
    0 0 0 0 0 1 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 7))
% assert(isequaln(obs.i_next, 8))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 3))
assert(isequaln(obs.f_main, [1 4 2 5]))

% Check probabilities
assert(isequal(obs.rk, [1 1 2 1 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.1865 0.6345 0.1865 0.8015 0.1865]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9638 0.0039 0.0097 0.0097 0.0097]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.9281 0.0129 0.0094 0.0403 0.0094]'))

% For comparison: probability densities if yk had been 0:
% assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [3.1646 1.2260 3.1646 3.1646 3.1646]'))
% assert(isequal(round(obs.p_seq_g_Yk, 4), [0.0050    0.0019    0.0050    0.4940    0.4940]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [0.4290; 0.1510]))
assert(isequaln(round(obs.xk_est, 4), [0.3971; 0.1510]))
assert(isequaln(round(obs.yk_est, 4), 0.1191))
assert(isequal(round(obs.yk_est - yk, 4), -0.1809))

% Update at k = 7
i = 8;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0
    0 0 0 0 0 0 1 0
    0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 8))
% assert(isequaln(obs.i_next, 1))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 5))
assert(isequaln(obs.f_main, [1 4 2 3]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 2]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.99 0.99 0.01]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0166 0.9494 0.0166 1.9388 0.0166]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.9188 0.0127 0.0093 0.0399 0.0093]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.1454 0.1150 0.0015 0.7366 0.0015]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [1.8619; 0.8147]))
assert(isequaln(round(obs.xk_est, 4), [1.4961; 0.8147]))
assert(isequaln(round(obs.yk_est, 4), 0.4488))
assert(isequal(round(obs.yk_est - yk, 4), -0.0612))

% Update at k = 8
i = 9;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
    0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0
    1 0 0 0 1 0 0 0  % New additions have looped back to 1st position
    0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 1))
% assert(isequaln(obs.i_next, 2))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, 3))
assert(isequaln(obs.f_main, [4 1 2 5]))  % most probable seq. has changed

% Check probabilities
assert(isequal(obs.rk, [1 1 2 1 1]'))
assert(isequal(obs.p_rk_g_rkm1, [0.99 0.99 0.01 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.0084 1.3806 2.5317 2.5317 0.0084]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.1439 0.1139 0.0074 0.7293 0.0015]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), [0.0006 0.0777 0.0092 0.9125 0]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [2.4556; 0.9601]))
assert(isequaln(round(obs.xk_est, 4), [2.1365; 0.9601]))
assert(isequaln(round(obs.yk_est, 4), 0.6409))
assert(isequal(round(obs.yk_est - yk, 4), -0.0161))


%% Run full simulation - SISO system

clear all
plot_dir = 'plots';

seed = 0;
rng(seed)

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Simulation settings
nT = 100;
t = Ts*(0:nT)';

% Choose time and amplitude of input disturbance
t_shock = 10;
du0 = 1;

% Inputs
%U = (idinput(size(t)) + 1)/2;
U = zeros(size(t));
U(t >= 1) = -1;
% Disturbance input
% This is used by the SKF observer
alpha = zeros(nT+1, 1);
alpha(t == t_shock(1), 1) = 1;
Wp = du0' .* alpha;

% Calculate the input disturbance
P = zeros(size(U));
P(t >= t_shock, 1) = du0;

% Combined inputs for simulation
U_sim = [U Wp];

% Simulate system
X = zeros(nT+1,n);
Y = zeros(nT+1,ny);
xk = zeros(n,1);

for i = 1:nT+1

    % Inputs
    uk = U_sim(i, :)';

    % Compute y(k)
    yk = C*xk + D*uk;

    % Store results
    X(i, :) = xk';
    Y(i, :) = yk';
    
    % Compute x(k+1)
    xk = A*xk + B*uk;

end

% Check simulation output is correct
[Y2, t, X2] = lsim(Gpss,U_sim,t);
assert(isequal(X, X2))
assert(isequal(Y, Y2))

% Choose measurement noise for plant
sigma_MP = 0;  % Set to zero for testing
Y_m = Y + sigma_MP'.*randn(size(Y));

% Define custom MKF test observers

% Devise a custom multi-model filter with a shock indicator 
% sequence that perfectly reflects the shock occurence in
% this test simulation (t = t_shock)
% Multiple model filter 1
obs_model = model;
obs_model.B = Bu;
obs_model.D = Du;
obs_model.R = sigma_M^2;
models = {obs_model, obs_model};
P0 = 1000*eye(n);
Q0 = diag([q1 0]);
P0_init = repmat({P0}, 1, 2);
models{1}.Q = diag([Q0(1,1) sigma_wp{1}(1)^2]);
models{2}.Q = diag([Q0(1,1) sigma_wp{1}(2)^2]);
seq = {ones(1, nT+1); ones(1, nT+1)};
seq{2}(t == 10) = 2;
p_rk = [1-epsilon epsilon]';
T = repmat(p_rk', 2, 1);
MKF3 = MKFObserverS(models,P0,seq,T,"MKF3");

seq = {ones(1, nT+1)};
seq{1}(t == 10) = 2;
p_rk = [1-epsilon epsilon]';
T = repmat(p_rk', 2, 1);
MKF4 = MKFObserverS(models,P0,seq,T,"MKF4");

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = SKFObserverS(models,P0,seq{1},"SKF");

% Simulate observers

% Choose observers to test
observers = {KF2, KF3, SKF, MKF3, MKF4, MKF_SP1, MKF_SP2};
% Note: KF1 is too slow to pass static error test here

% Measured inputs (not including disturbances)
U_m = U;
n_obs = numel(observers);
MSE = struct();

for i = 1:n_obs

    obs = observers{i};
    [obs, sim_results] = run_test_simulation(nT,Ts,n,ny,U_m,Y_m,obs);

    % Check observer errors are zero prior to
    % input disturbance
    assert(all(abs(sim_results.X_est(1:20,:) - X(1:20, :)) < 1e-10, [1 2]))
    assert(all(abs(sim_results.Y_est(1:20,:) - Y(1:20, :)) < 1e-10))

    % Check observer static errors are small
    % after input disturbance
    if all(sigma_MP == 0)
        assert(abs(sim_results.Y_est(end, :) - Y(end, :)) < 3e-4);
        assert(abs(sim_results.X_est(end, 2) - du0) < 4e-4);
    end
    % TODO: Errors for AFMM1 were not as low as for RODD MKF observers

    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE.(obs.label) = mean((Y_est - Y).^2);

    % Save updated observer
    observers{i} = obs;

end

MSE_test_values = struct(...
    'MKF_SP1', 0.000491, ...
    'MKF_SP2', 0.000492, ...
    'KF2', 0.000006, ...
    'KF3', 0.000992, ...
    'SKF', 0.000012, ...
    'MKF3', 0.000501, ...
    'MKF4', 0.000012 ...
);
labels = fieldnames(MSE);
% for i = 1:numel(labels)
%     fprintf("%s: %f (%f)\n", labels{i}, MSE.(labels{i}), MSE_test_values.(labels{i}))
% end
for i = 1:numel(labels)
    assert(isequal(round(MSE.(labels{i}), 6), MSE_test_values.(labels{i})), ...
        labels{i})
end

% % Display results of last simulation
% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% table(t,alpha,U,P,Wp,X,Y,Y_m,X_est,Y_est,E_obs)
% 
% K_obs = sim_results.K_obs;
% trP_obs = sim_results.trP_obs;
% 
% % Display gains and trace of covariance matrix
% table(t, cell2mat(K_obs), trP_obs, ...
%     'VariableNames',{'t', 'K{1}, K{2}', 'trace(P)'})
% 
% % Display MKF_SP filter groupings
% switch obs.type
%     case {'MKF_SP'}
%     f_hold = sim_results.MKF_SP_f_hold
%     f_main = sim_results.MKF_SP_f_main
%     [array2table(f_hold) array2table(f_main)]
% end
% 
% % Plot of inputs and outputs
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
% ylabel('$y_m(k)$ and $\hat{y}(k)$')
% title('Process output measurements and estimates')
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
% ylabel('$x_i(k)$ and $\hat{x}_i(k)$')
% labels = repmat({''}, 1, n*2);
% for i = 1:n
%     labels{i} = sprintf("$x_{%d}(k)$", i);
% end
% for i = 1:n
%     labels{i+n} = sprintf("$%s{x}_{%d}(k)$", '\hat', i);
% end
% legend(labels)
% title('Actual states and estimates')
% grid on
% 
% ax3 = subplot(4,1,3);
% stairs(t,U,'Linewidth',2); hold on;
% stairs(t,Wp,'Linewidth',2)
% stairs(t,P,'Linewidth',2)
% max_min = [min(min([U Wp P])) max(max([U Wp P]))];
% bd = max([0.1 diff(max_min)*0.1]);
% ylim(max_min + [-bd bd])
% ylabel('$u(k)$, $w_p(k)$ and $p(k)$')
% legend('$u(k)$', '$w_p(k)$', '$p(k)$')
% title('Actual process inputs')
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
% 
% 
% % Plot of conditional filter probabilities
% switch obs.type
%     case {'MKF', 'MKF_S', 'MKF_SF', 'MKF_SP'}
%         p_seq_g_Yk = sim_results.MKF_p_seq_g_Yk;
%         % Note: first data points are nans,
%         % ignore last data point to make plot wider
% 
%         figure(11); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Pr(\Gamma(k) \mid Y(k))$'};
%         make_waterfall_plot(t(2:end-1), p_seq_g_Yk(2:end-1, :), [0 1], ...
%             ax_labels, [0 82]);
%         filename = sprintf('rod_mkf_observer_test_pyk_wfplot.pdf');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Conditional probabilities of y(k)')
% end
% 
% % Plot of trace of filter covariance matrices
% switch obs.type
%     case {'MKF', 'MKF_S', 'MKF_SP'}
%         figure(12); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Tr(P(k))$'};
%         make_waterfall_plot(t, sim_results.trP_obs_f, [0 5], ax_labels, [0 82]);
%         filename = sprintf('rod_mkf_observer_test_trP_wfplot.pdf');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Trace of covariance matrices')
% end
% 
% % Plot of filter correction gains (k1)
% switch obs.type
%     case {'MKF', 'MKF_S', 'MKF_SP'}
%         K_obs = cell2mat(sim_results.K_obs_f);
%         % Select first gain value onlu
%         K1_obs = K_obs(:,1:2:end);
% 
%         figure(13); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$K_1$'};
%         make_waterfall_plot(t, K1_obs, [0 0.5], ax_labels, [0 82]);
%         filename = sprintf('rod_mkf_observer_test_K_wfplot.png');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Filter correction gains (k1)')
% end
% 
% % Plot of final sequence values
% switch obs.type
%     case 'MKF_SPS'
%         Z = double(cell2mat(obs.seq))';
%         if size(Z, 1) > nT
%             Z = Z(1:nT,:);
%         else
%             Z = [Z(1:obs.i,:); Z(obs.i+1:end,:)];
%         end
%         seq_len = size(Z, 1);
%         t = Ts*(nT-seq_len+1:nT)';
% 
%         figure(14); clf
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$\gamma(k)$'};
%         title('Final filter sequence values')
%         make_waterfall_plot(t,Z,[0 1], ax_labels, [0 82]);
%         filename = sprintf('rod_afmm_filter_test.png');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
% end


%% Test initialization on 2x2 system

% Load system and disturbance model from file
sys_rodin_step_2x2sym

% Load observers from file
obs_rodin_step_2x2

% Check observer initialization
assert(isequal(MKF_SP1.epsilon, epsilon))
assert(isequal(MKF_SP1.sigma_wp, sigma_wp))
assert(isequal(MKF_SP1.Q0, Q0))
assert(isequal(MKF_SP1.R, R))
assert(MKF_SP1.nh == 19)
assert(MKF_SP1.n_min == 5)
assert(isequal(MKF_SP1.n_hold, 5*2))
assert(isequal(MKF_SP1.n_main, 9))
assert(isequaln(MKF_SP1.f_hold, 10:19))
assert(isequaln(MKF_SP1.f_main, 1:9))
% assert(isequaln(MKF_SP1.i, 0))
assert(MKF_SP1.n == 4)
assert(MKF_SP1.nu == 2)
assert(MKF_SP1.ny == 2)
assert(MKF_SP1.nj == 3)
assert(isequal(MKF_SP1.sys_model, model))
assert(MKF_SP1.Ts == Ts)
assert(isequaln(MKF_SP1.io.u_known, u_known))
assert(isequal(MKF_SP1.models{1}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(1)^2])))
assert(isequal(MKF_SP1.models{2}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(2)^2 sigma_wp{2}(1)^2])))
assert(isequal(MKF_SP1.models{3}.Q, ...
    diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(2)^2])))
assert(isequal([MKF_SP1.models{1}.R MKF_SP1.models{2}.R MKF_SP1.models{3}.R], ...
    repmat(R, 1, 3)))
%assert(isequal(size(MKF_SP1.seq), [MKF_SP1.nh 1]))
%assert(isequal(size(cell2mat(MKF_SP1.seq)), [MKF_SP1.nh MKF_SP1.f]))
%assert(MKF_SP1.f == size(MKF_SP1.seq{1}, 2))
assert(isequal(MKF_SP1.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SP1.Pkp1, 1000*eye(4)))
assert(isequal(MKF_SP1.r0, ones(MKF_SP1.nh, 1)))
assert(isequal(MKF_SP1.p_seq_g_Yk_init, [1; zeros(MKF_SP1.nh-1, 1)]))
assert(isequaln(MKF_SP1.p_rk_g_rkm1, nan(MKF_SP1.nh, 1)))
assert(isequaln(MKF_SP1.xk_est, nan(n, 1)))
assert(isequaln(MKF_SP1.Pk, nan(n, n)))
assert(isequaln(MKF_SP1.yk_est, nan(ny, 1)))

% Check observer initialization
assert(isequal(MKF_SP2.epsilon, epsilon))
assert(isequal(MKF_SP2.sigma_wp, sigma_wp))
assert(isequal(MKF_SP2.Q0, Q0))
assert(isequal(MKF_SP2.R, R))
assert(MKF_SP2.nh == 25)
assert(MKF_SP2.n_min == 9)
assert(isequal(MKF_SP2.n_hold, 9*2))
assert(isequal(MKF_SP2.n_main, 7))
assert(isequaln(MKF_SP2.f_hold, 8:25))
assert(isequaln(MKF_SP2.f_main, 1:7))
% assert(isequaln(MKF_SP2.i, 0))
assert(MKF_SP2.n == 4)
assert(MKF_SP2.nu == 2)
assert(MKF_SP2.ny == 2)
assert(MKF_SP2.nj == 3)
assert(isequal(MKF_SP1.sys_model, model))
assert(MKF_SP2.Ts == Ts)
assert(isequaln(MKF_SP2.io.u_known, u_known))
assert(isequal(MKF_SP2.models{1}.Q, diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(1)^2])))
assert(isequal(MKF_SP2.models{2}.Q, diag([0.01 0.01 sigma_wp{1}(2)^2 sigma_wp{2}(1)^2])))
assert(isequal(MKF_SP2.models{3}.Q, diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(2)^2])))
assert(isequal([MKF_SP2.models{1}.R MKF_SP1.models{2}.R MKF_SP1.models{3}.R], ...
    repmat(R, 1, 3)))
% assert(isequal(size(MKF_SP2.seq), [MKF_SP2.nh 1]))
% assert(isequal(size(cell2mat(MKF_SP2.seq)), [MKF_SP2.nh MKF_SP2.f]))
% assert(MKF_SP2.f == size(MKF_SP2.seq{1}, 2))
assert(isequal(MKF_SP2.xkp1_est, zeros(n, 1)))
assert(isequal(MKF_SP2.Pkp1, 1000*eye(4)))
assert(isequal(MKF_SP2.r0, ones(MKF_SP2.nh, 1)))
assert(isequal(MKF_SP2.p_seq_g_Yk_init, [1; zeros(MKF_SP2.nh-1, 1)]))
assert(isequaln(MKF_SP2.p_rk_g_rkm1, nan(MKF_SP2.nh, 1)))
assert(isequaln(MKF_SP2.xk_est, nan(n, 1)))
assert(isequaln(MKF_SP2.Pk, nan(n, n)))
assert(isequaln(MKF_SP2.yk_est, nan(ny, 1)))

% Check optional definition with an initial state estimate works
x0 = [0.1; 0.5; -0.2; -0.4];
io.u_known = u_known;
io.y_meas = true(ny, 1);
d = 1;
MKF_SP_testx0 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,MKF_SP2.nh,n_min,label,d,x0);
assert(isequal(MKF_SP_testx0.xkp1_est, x0))
assert(isequal(MKF_SP_testx0.r0, ones(MKF_SP_testx0.nh, 1)))
assert(isequal(MKF_SP_testx0.p_seq_g_Yk_init, [1; zeros(MKF_SP_testx0.nh-1, 1)]))


%% Test sequence on 2x2 system

% Load system and disturbance model from file
sys_rodin_step_2x2sym

% Values use in this test may be different to above
epsilon = [0.005 0.005]';
sigma_wp = {[0.0100 1], [0.0100 1]};
sigma_M = [0.1 0.1]';

% Load observers from file
obs_rodin_step_2x2

x0 = [0; 0; 0; 0];
nh = 10;  % number of filters
n_min = 3;  % minimum life of cloned filters
io.u_known = u_known;
io.y_meas = true(ny, 1);
d = 1;
obs = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nh,n_min, ...
    label,d,x0);
assert(isequal(obs.xkp1_est, x0))

% Test simuation data
sim_data = [ ...
         0         0         0         0         0         0         0;
    1.0000         0         0         0         0         0         0;
    2.0000         0         0         0         0         0         0;
    3.0000         0         0         0         0         0         0;
    4.0000         0         0    1.0000         0         0         0;
    5.0000         0         0         0         0         0         0;
    6.0000         0         0         0         0    0.1110   -0.0222;
    7.0000         0         0    1.0000         0    0.2097   -0.0419;
    8.0000         0         0         0         0    0.2974   -0.0595;
    9.0000         0         0         0         0    0.4864   -0.0973;
   10.0000         0         0         0         0    0.6544   -0.1309
];
nT = size(sim_data, 1);
t = sim_data(:, 1);
U_m = sim_data(:, 2:3);
Wp = sim_data(:, 4:5);
Y_m = sim_data(:, 6:7);

% Set marker values on each sequence - for testing only
% these values at the end of the sequences are not used
% by the observer.
% for i = 1:nh
%     obs.seq{i}(8) = i;
% end
seq0 = [
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   2
   0   0   0   0   0   0   0   3
   0   0   0   0   0   0   0   4
   0   0   0   0   0   0   0   5
   0   0   0   0   0   0   0   6
   0   0   0   0   0   0   0   7
   0   0   0   0   0   0   0   8
   0   0   0   0   0   0   0   9
   0   0   0   0   0   0   0  10
];
%%disp(obs.i)  % use for debugging
%disp(debug_array(obs))
% assert(isequaln(obs.i, 0))
% assert(isequaln(obs.i_next, 1))
% assert(isequaln(cell2mat(obs.seq), seq0))
assert(isequal(obs.n_hold, 6))
assert(isequal(obs.n_main, 4))
assert(isequaln(obs.f_hold, [5 6 7 8 9 10]))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 1 1 1 1 1 1]'))
assert(isequaln(round(obs.p_rk_g_rkm1, 4), nan(10,1)))
assert(isequaln(obs.p_yk_g_seq_Ykm1, nan(10,1)))
assert(isequaln(obs.p_seq_g_Ykm1, nan(10,1)))
assert(isequal(obs.p_seq_g_Yk, [1 zeros(1, 9)]'))

% Check initialization of filters
assert(isequal(obs.filters.Pkp1(:,:,1), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,2), 1e10*eye(4)))
assert(isequal(obs.filters.Pkp1(:,:,3), 1e10*eye(4)))
assert(isequal(obs.filters.Pkp1(:,:,4), 1e10*eye(4)))
assert(isequal(obs.filters.Pkp1(:,:,5), 1e10*eye(4)))

% Check estimates
assert(isequal(obs.xkp1_est, [0 0 0 0]'))
assert(isequaln(obs.xk_est, nan(4, 1)))
assert(isequaln(obs.yk_est, nan(2, 1)))

% Update at k = 0
i = 1;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   2
   1   0   0   0   0   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   5
   0   0   0   0   0   0   0   6
   0   0   0   0   0   0   0   7
   0   0   0   0   0   0   0   8
   0   0   0   0   0   0   0   9
   0   0   0   0   0   0   0  10
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 1))
% assert(isequaln(obs.i_next, 2))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [7 8 9 10 3 4]))
assert(isequaln(obs.f_main, [1 2 5 6]))

% Check probabilities
assert(isequal(obs.rk, [1 1 2 3 1 1 1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.99 0.005 0.005 0.99 0.99 0.99 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0129 0 0.0129 0.0129 0 0 0 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.99 0 0.005 0.005 0 0 0 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.99 0 0.005 0.005 0 0 0 0 0 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0 0 0 0]'))
assert(isequaln(obs.xk_est, [0 0 0 0]'))
assert(isequaln(obs.yk_est, [0 0]'))

% Update at k = 1
i = 2;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   2
   1   0   0   0   0   0   0   1
   2   0   0   0   0   0   0   1
   0   1   0   0   0   0   0   1
   0   2   0   0   0   0   0   1
   0   0   0   0   0   0   0   7
   0   0   0   0   0   0   0   8
   0   0   0   0   0   0   0   9
   0   0   0   0   0   0   0  10
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 2))
% assert(isequaln(obs.i_next, 3))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [9 10 3 4 5 6]))
assert(isequaln(obs.f_main, [1 2 7 8]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 2 3 1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.99 0.99 0.99 0.005 0.005 0.99 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0172 0 0.0172 0.0172 0.0172 0.0172 0 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9802 0 0.0049 0.0049 0.0049 0.0049 0 0 0 0]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9803 0 0.0049 0.0049 0.0049 0.0049 0 0 0 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0 0 0 0]'))
assert(isequaln(obs.xk_est, [0 0 0 0]'))
assert(isequaln(obs.yk_est, [0 0]'))

% Update at k = 2
i = 3;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   2
   1   0   0   0   0   0   0   1
   2   0   0   0   0   0   0   1
   0   1   0   0   0   0   0   1
   0   2   0   0   0   0   0   1
   0   0   1   0   0   0   0   1
   0   0   2   0   0   0   0   1
   0   0   0   0   0   0   0   9
   0   0   0   0   0   0   0  10
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 3))
% assert(isequaln(obs.i_next, 4))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [3 4 5 6 7 8]))
assert(isequaln(obs.f_main, [1 2 9 10]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 1 1 2 3 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.99 0.99 0.99 0.99 0.99 0.005 0.005 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [2.9639 2.9564 2.6133 2.6133 2.9639 2.9639 2.9639 2.9639 2.9564 2.9564]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9705 0 0.0049 0.0049 0.0049 0.0049 0.0049 0.0049 0 0]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9719 0 0.0043 0.0043 0.0049 0.0049 0.0049 0.0049 0 0]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0 0 0 0]'))
assert(isequaln(obs.xk_est, [0 0 0 0]'))
assert(isequaln(obs.yk_est, [0 0]'))

% Update at k = 3
i = 4;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   2
   1   0   0   0   0   0   0   1
   2   0   0   0   0   0   0   1
   0   1   0   0   0   0   0   1
   0   2   0   0   0   0   0   1
   0   0   1   0   0   0   0   1
   0   0   2   0   0   0   0   1
   0   0   0   1   0   0   0   1
   0   0   0   2   0   0   0   1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 4))
% assert(isequaln(obs.i_next, 5))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [5 6 7 8 9 10]))
assert(isequaln(obs.f_main, [1 2 3 4]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 1 1 1 1 2 3]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.005 0.005]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [5.3049 5.3019 4.8905 4.8905 4.3126 4.3126 5.3049 5.3049 5.3049 5.3049]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9622 0 0.0043 0.0043 0.0048 0.0048 0.0048 0.0048 0.0048 0.0048]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9648 0 0.0039 0.0039 0.0039 0.0039 0.0048 0.0048 0.0048 0.0048]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0 0 0 0]'))
assert(isequaln(obs.xk_est, [0 0 0 0]'))
assert(isequaln(obs.yk_est, [0 0]'))

% Update at k = 4
i = 5;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   2   0   0   1  % these seq.s have now been replaced
   0   0   0   0   1   0   0   1  %
   2   0   0   0   0   0   0   1
   0   1   0   0   0   0   0   1
   0   2   0   0   0   0   0   1
   0   0   1   0   0   0   0   1
   0   0   2   0   0   0   0   1
   0   0   0   1   0   0   0   1
   0   0   0   2   0   0   0   1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 5))
% assert(isequaln(obs.i_next, 6))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [7 8 9 10 3 2]))
assert(isequaln(obs.f_main, [1 4 5 6]))

% Check probabilities
assert(isequal(obs.rk, [1 3 2 1 1 1 1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.005 0.005 0.99 0.99 0.99 0.99 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [7.0447 7.0447 7.0447 6.6822 5.8568 5.8568 5.4327 5.4327 7.0447 7.0447]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9552 0.0048 0.0048 0.0039 0.0039 0.0039 0.0048 0.0048 0.0048 0.0048]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9629 0.0048 0.0048 0.0037 0.0033 0.0033 0.0037 0.0037 0.0048 0.0048]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0 0 0 0]'))
assert(isequaln(obs.xk_est, [0 0 0 0]'))
assert(isequaln(obs.yk_est, [0 0]'))

% Update at k = 5
i = 6;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   2   0   0   1
   0   0   0   0   1   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   0   2   0   1
   0   0   0   0   0   1   0   1
   0   0   1   0   0   0   0   1
   0   0   2   0   0   0   0   1
   0   0   0   1   0   0   0   1
   0   0   0   2   0   0   0   1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 6))
% assert(isequaln(obs.i_next, 7))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [9 10 3 2 6 5]))
assert(isequaln(obs.f_main, [1 4 7 8]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 3 2 1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.99 0.99 0.99 0.005 0.005 0.99 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [8.3545 8.3545 8.3545 8.0562 8.3545 8.3545 6.4384 6.4384 6.213 6.213]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9533 0.0048 0.0048 0.0037 0.0048 0.0048 0.0037 0.0037 0.0048 0.0048]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9641 0.0048 0.0048 0.0036 0.0048 0.0048 0.0029 0.0029 0.0036 0.0036]'))

% Check estimates
assert(isequal(obs.xkp1_est, [0 0 0 0]'))
assert(isequaln(obs.xk_est, [0 0 0 0]'))
assert(isequaln(obs.yk_est, [0 0]'))

% Update at k = 6  *** First non-zero measurement ***
i = 7;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   1
   0   0   0   0   2   0   0   1
   0   0   0   0   1   0   0   1
   2   0   0   0   0   0   0   1
   0   0   0   0   0   2   0   1
   0   0   0   0   0   1   0   1
   0   0   0   0   0   0   2   1
   0   0   0   0   0   0   1   1
   0   0   0   1   0   0   0   1
   0   0   0   2   0   0   0   1
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 7))
% assert(isequaln(obs.i_next, 8))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [3 2 6 5 8 7]))
assert(isequaln(obs.f_main, [1 4 9 10]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 1 1 3 2 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.99 0.99 0.99 0.99 0.99 0.005 0.005 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [6.4234 4.9777 5.4971 6.3022 6.4234 6.4234 6.4234 6.4234 5.5241 5.0106]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9544 0.0048 0.0048 0.0036 0.0048 0.0048 0.0048 0.0048 0.0036 0.0036]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9633 0.0038 0.0041 0.0035 0.0048 0.0048 0.0048 0.0048 0.0031 0.0028]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [0.4793 -0.0994 0.1322 0.0418]'))
assert(isequaln(round(obs.xk_est, 4), [0.4140 -0.0845 0.1322 0.0418]'))
assert(isequaln(round(obs.yk_est, 4), [0.0460 -0.0094]'))
assert(isequal(round(obs.yk_est - yk, 4), [-0.0650 0.0128]'))

% Update at k = 7
i = 8;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   0
   0   0   0   0   2   0   0   0
   0   0   0   0   1   0   0   0
   2   0   0   0   0   0   0   0
   0   0   0   0   0   2   0   0
   0   0   0   0   0   1   0   0
   0   0   0   0   0   0   2   0
   0   0   0   0   0   0   1   0
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   2
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 8))
% assert(isequaln(obs.i_next, 1))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [6 5 8 7 9 10]))
assert(isequaln(obs.f_main, [1 4 3 2]))

% Check probabilities
assert(isequal(obs.rk, [1 1 1 1 1 1 1 1 2 3]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.005 0.005]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [4.4821 4.1027 6.0731 4.5274 3.7076 4.6422 4.4821 4.4821 4.4821 4.4821]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9537 0.0037 0.0041 0.0035 0.0048 0.0048 0.0048 0.0048 0.0048 0.0048]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9592 0.0034 0.0056 0.0036 0.004 0.005 0.0048 0.0048 0.0048 0.0048]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [1.1179 -0.2320 0.2794 0.0881]'))
assert(isequaln(round(obs.xk_est, 4), [0.9928 -0.2029 0.2794 0.0881]'))
assert(isequaln(round(obs.yk_est, 4), [0.1102 -0.0225]'))
assert(isequal(round(obs.yk_est - yk, 4), [-0.0995 0.0194]'))

% Update at k = 8
i = 9;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
seq = [
   0   0   0   0   0   0   0   0
   2   0   0   0   0   0   0   0  % New additions have looped back to 1st position
   0   0   0   0   1   0   0   0
   1   0   0   0   0   0   0   0  %
   0   0   0   0   0   2   0   0
   0   0   0   0   0   1   0   0
   0   0   0   0   0   0   2   0
   0   0   0   0   0   0   1   0
   0   0   0   0   0   0   0   1
   0   0   0   0   0   0   0   2
];
%disp(obs.i)
%display_obs_data(obs)
%disp(debug_array(obs))
% assert(isequaln(obs.i, 1))
% assert(isequaln(obs.i_next, 2))
% assert(isequaln(cell2mat(obs.seq), seq))
assert(isequaln(obs.f_hold, [8 7 9 10 4 2]))
assert(isequaln(obs.f_main, [1 3 6 5]))

% Check probabilities
assert(isequal(obs.rk, [1 3 1 2 1 1 1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.99 0.005 0.99 0.005 0.99 0.99 0.99 0.99 0.99 0.99]'))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [3.6803 3.6803 7.4531 3.6803 3.6741 6.465 3.1649 4.2927 3.6803 3.6803]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.9496 0.0048 0.0055 0.0048 0.0039 0.0049 0.0048 0.0048 0.0048 0.0048]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9475 0.0048 0.0112 0.0048 0.0039 0.0087 0.0041 0.0056 0.0048 0.0048]'))

% Check estimates
assert(isequal(round(obs.xkp1_est, 4), [1.8033 -0.3779  0.4175  0.1298]'))
assert(isequaln(round(obs.xk_est, 4), [1.6319 -0.3363 0.4175 0.1298]'))
assert(isequaln(round(obs.yk_est, 4), [0.1811 -0.0373]'))
assert(isequal(round(obs.yk_est - yk, 4), [-0.1163 0.0222]'))


%% Full simulation on 2x2 system

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
Gpss = ss(A,B,C,D,Ts);

% Dimensions
n = size(A, 1);
ny = size(C, 1);

% Model parameter struct used by observers
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.D = D;
model.Ts = Ts;

% Designate measured input and output signals
u_known = [true; true; false; false];
y_meas = [true; true];

% Observer model without disturbance noise input
Bu = B(:, u_known);
Du = D(:, u_known);
nu = sum(u_known);
nw = sum(~u_known);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_known);
nw = sum(~u_known);

% RODD random variable parameters
epsilon = [0.01; 0.01];
sigma_M = [0.1; 0.1];
sigma_wp = {[0.01 1], [0.01 1]};

% Different values for covariance matrix
Q1 = diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(1)^2]);
Q2 = diag([0.01 0.01 sigma_wp{1}(2)^2 sigma_wp{2}(1)^2]);
Q3 = diag([0.01 0.01 sigma_wp{1}(1)^2 sigma_wp{2}(2)^2]);

% Covariance of output errors
R = diag(sigma_M.^2);

% Observer models for new observer functions
models = {struct, struct};
models{1}.A = A;
models{1}.B = Bu;
models{1}.C = C;
models{1}.Ts = Ts;
models{1}.Q = Q1;
models{1}.R = R;
models{2}.A = A;
models{2}.B = Bu;
models{2}.C = C;
models{2}.Ts = Ts;
models{2}.Q = Q2;
models{2}.R = R;
models{3}.A = A;
models{3}.B = Bu;
models{3}.C = C;
models{3}.Ts = Ts;
models{3}.Q = Q3;
models{3}.R = R;

% Kalman filter 3 - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
model_KF3 = models{3};
model_KF3.Q = diag([0.01 0.01 0.1^2 0.1^2]);
KF3 = KalmanFilterF(model_KF3,P0,'KF3');

% Multiple model observer with sequence pruning 1
label = "MKF_SP1";
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
nh = 15;  % number of filters
n_min = 5;  % minimum life of cloned filters
io.u_known = u_known;
io.y_meas = true(ny, 1);
MKF_SP1 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nh, ...
    n_min,label);

% Multiple model observer with sequence pruning 2
label = "MKF_SP2";
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
nh = 30;  % number of filters
n_min = 10;  % minimum life of cloned filters
io.u_known = u_known;
io.y_meas = true(ny, 1);
MKF_SP2 = MKFObserverSP_RODD(model,io,P0,epsilon,sigma_wp, ...
    Q0,R,nh,n_min,label);

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
A2 = repmat({A}, 1, 3);
Bu2 = repmat({Bu}, 1, 3);
C2 = repmat({C}, 1, 3);
Du2 = repmat({Du}, 1, 3);
P0 = 1000*eye(n);
%P0_init = repmat({P0}, 1, 3);

seq = repmat({ones(1, nT+1)}, 4, 1);
seq{2}(t == t_shock(1)) = 2;  % shock 1
seq{3}(t == t_shock(2)) = 3;  % shock 2
seq{4}(t == t_shock(1)) = 2;  % both
seq{4}(t == t_shock(2)) = 3;
p_rk = [1-epsilon epsilon]';
Z = [1 1; 2 1; 1 2];  % combinations
p_rk = prod(prob_rk(Z', p_rk), 1)';
p_rk = p_rk ./ sum(p_rk);  % normalized
T = repmat(p_rk', 3, 1);
MKF3 = MKFObserverS(models,P0,seq,T,'MKF3');
assert(MKF3.nh == 4)

seq = {ones(1, nT+1)};
seq{1}(t == t_shock(1)) = 2;
seq{1}(t == t_shock(2)) = 3;
MKF4 = MKFObserverS(models,P0,seq,T,'MKF4');
assert(MKF4.nh == 1)

% Define scheduled Kalman filter
% Note: in the case of more than one random input variable, all
% possible combinations of the switching systems need to be 
% accounted for.
% Here, we account for 3 possible combinations:
% combs = [0 0; 1 0; 0 1];
% (This is the same as the MKF filters for the RODD).
% seq = sum(alpha .* 2.^(1:-1:0), 2)';
SKF = SKFObserverS(models,P0,seq{1},"SKF");

% Choose observers to test
observers = {KF3, MKF_SP1, MKF_SP2, MKF3, MKF4, SKF};

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
MSE = struct();
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
    MSE.(obs.label) = mean((Y_est - Y).^2);
    %fprintf("%d, %s: %f\n", i, obs.label, mean((Y_est - Y).^2))

    % Save updated observer
    observers{i} = obs;

end


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

% Results on Nov 8 after reverting back the Bayesian updating
MSE_test_values = struct(...
    'KF3', [0.000296 0.000433], ...
    'MKF_SP1', [0.000308 0.000758], ...
    'MKF_SP2', [0.000322 0.000773], ...
    'MKF3', [0.000332 0.000338], ...
    'MKF4', [0.000017 0.000022], ...
    'SKF', [0.000017 0.000022] ...
);

labels = fieldnames(MSE);
% for i = 1:numel(labels)
%     fprintf("%s: %f %f (%f %f)\n", labels{i}, MSE.(labels{i}), ...
%         MSE_test_values.(labels{i}))
% end
for i = 1:numel(labels)
    assert(isequal(round(MSE.(labels{i}), 6), MSE_test_values.(labels{i})), ...
        labels{i})
end

% END


% % Functions for debugging and testing sequences
% function display_obs_data(obs)
%     disp(obs.f_hold); disp(obs.f_main); 
%     disp(strjoin(string(round(obs.rk, 4)')));
%     disp(strjoin(string(round(obs.p_rk_g_rkm1, 4)')));
%     disp(strjoin(string(round(obs.p_yk_g_seq_Ykm1, 4)')));
%     disp(strjoin(string(round(obs.p_seq_g_Ykm1, 4)')));
%     disp(strjoin(string(round(obs.p_seq_g_Yk, 4)')));
% end
% 
% function dba = debug_array(obs)
%     hold = zeros(obs.nh, 1);
%     hold(nonzeros(obs.f_hold)) = nonzeros(obs.f_hold);
%     main = zeros(obs.nh, 1);
%     main(nonzeros(obs.f_main)) = nonzeros(obs.f_main);
%     seq = cell2mat(obs.seq);
%     p_max = (obs.p_seq_g_Yk == max(obs.p_seq_g_Yk));
%     dba = [table(hold, main) array2table(seq) table(p_max)];
% end