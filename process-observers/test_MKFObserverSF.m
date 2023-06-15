% Tests the following classes:
% 
%  - MKFObserverSF
%  - MKFObserverSF_DI
%  - MKFObserverSF_RODD
%

clear all

% If plotting:
plot_dir = 'plots';
addpath("../plot-utils")

seed = 0;
rng(seed)


%% Test initialization with SISO system

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Check observer attributes
assert(strcmp(MKF_SF1.type, "MKF_SF_RODD"))
assert(strcmp(MKF_SF1.label, "MKF_SF1"))
assert(isequal(MKF_SF1.sys_model, model))
assert(isequal(MKF_SF1.io, io))
assert(isequal(MKF_SF1.P0, P0))
assert(isequal(MKF_SF1.P0, P0))
assert(MKF_SF1.epsilon == epsilon)
assert(isequal(MKF_SF1.sigma_wp, sigma_wp))
assert(isequal(MKF_SF1.Q0, Q0))
assert(isequal(MKF_SF1.R, R))
assert(isequal(MKF_SF1.nf, 3))
assert(isequal(MKF_SF1.m, 1))
assert(isequal(MKF_SF1.d, 5))
seq = { ...
    [1 1 1]
    [2 1 1]
    [1 2 1]
    [1 1 2]
};
assert(isequal(MKF_SF1.seq, seq))
T = [0.9510    0.0490;
     0.9510    0.0490];
assert(isequal(round(MKF_SF1.T, 4), T))
assert(MKF_SF1.nh == 6)
assert(MKF_SF1.nm == 4)
assert(MKF_SF1.nh_max == 6)
assert(size(MKF_SF1.filters.Xkp1_est, 3) == MKF_SF1.nh)
assert(isequaln(MKF_SF1.i, 3))  % last position in seq
assert(isequal(MKF_SF1.i_next, 1))
assert(MKF_SF1.n == n)
assert(MKF_SF1.nu == nu)
assert(MKF_SF1.ny == ny)
assert(MKF_SF1.nf == size(MKF_SF1.seq{1}, 2))
assert(MKF_SF1.nj == 2)
%assert(isequal(MKF_SF1.xkp1_est, zeros(n, 1)))
%assert(isequal(MKF_SF1.Pkp1, P0))
%assert(isequal(MKF_SF1.r0, [1 1 1 1 2 2]'))
assert(isequal(MKF_SF1.r0, [1 1 1 1 1 2]'))
assert(isequal(MKF_SF1.p_seq_g_Yk_init, ones(6, 1) ./ 6))
assert(isequal(MKF_SF1.rk, MKF_SF1.r0))
assert(isequal(MKF_SF1.rkm1, zeros(6, 1)))
assert(isequaln(MKF_SF1.p_yk_g_seq_Ykm1, nan(6, 1)))
assert(isequaln(MKF_SF1.p_rk_g_Ykm1, nan(6, 1)))
assert(isequaln(MKF_SF1.p_rk_g_rkm1, nan(6, 1)))
assert(isequaln(MKF_SF1.p_seq_g_Ykm1, nan(6, 1)))
assert(isequaln(MKF_SF1.xk_est, nan(2, 1)))
assert(isequaln(MKF_SF1.yk_est, nan))

% Check observer attributes
assert(strcmp(MKF_SF2.type, "MKF_SF_RODD"))
assert(strcmp(MKF_SF2.label, "MKF_SF2"))
assert(isequal(MKF_SF2.sys_model, model))
assert(isequal(MKF_SF2.io, io))
assert(isequal(MKF_SF2.P0, P0))
assert(isequal(MKF_SF2.P0, P0))
assert(MKF_SF2.epsilon == epsilon)
assert(isequal(MKF_SF2.sigma_wp, sigma_wp))
assert(isequal(MKF_SF2.Q0, Q0))
assert(isequal(MKF_SF2.R, R))
assert(isequal(MKF_SF2.nf, 5))
assert(isequal(MKF_SF2.m, 2))
assert(isequal(MKF_SF2.d, 3))
seq = { ...
    [1 1 1 1 1]
    [2 1 1 1 1]
    [1 2 1 1 1]
    [1 1 2 1 1]
    [1 1 1 2 1]
    [1 1 1 1 2]
    [2 2 1 1 1]
    [2 1 2 1 1]
    [2 1 1 2 1]
    [2 1 1 1 2]
    [1 2 2 1 1]
    [1 2 1 2 1]
    [1 2 1 1 2]
    [1 1 2 2 1]
    [1 1 2 1 2]
    [1 1 1 2 2]
};
assert(isequal(MKF_SF2.seq, seq))
T = [0.9703 0.0297
     0.9703 0.0297];
assert(isequal(round(MKF_SF2.T, 4), T))
assert(MKF_SF2.nh == 26)
assert(MKF_SF2.nm == 16)
assert(MKF_SF2.nh_max == 26)
assert(size(MKF_SF2.filters.Xkp1_est, 3) == MKF_SF2.nh)
assert(isequaln(MKF_SF2.i, 5))  % last position in seq
assert(isequal(MKF_SF2.i_next, 1))
assert(MKF_SF2.n == n)
assert(MKF_SF2.nu == nu)
assert(MKF_SF2.ny == ny)
assert(MKF_SF2.nf == size(MKF_SF2.seq{1}, 2))
assert(MKF_SF2.nj == 2)
%assert(isequal(MKF_SF2.xkp1_est, zeros(n, 1)))
%assert(isequal(MKF_SF2.Pkp1, P0))
% assert(isequal(MKF_SF2.r0, [ ...
%     1 1 1 1 1 1 1 1 1 1 ...
%     2 2 1 1 1 2 2 1 1 2 ...
%     2 1 2 2 2 2]' ...
% ))
assert(isequal(MKF_SF2.r0, [ ...
    1 1 1 1 1 1 1 1 1 1 ...
    2 2 1 1 1 1 1 1 2 2 ...
    1 1 2 1 2 2]' ...
))
assert(isequal(MKF_SF2.p_seq_g_Yk_init, ones(26, 1) ./ 26))
assert(isequal(MKF_SF2.rk, MKF_SF2.r0))
assert(isequal(MKF_SF2.rkm1, zeros(26, 1)))
assert(isequaln(MKF_SF2.p_yk_g_seq_Ykm1, nan(26, 1)))
assert(isequaln(MKF_SF2.p_rk_g_Ykm1, nan(26, 1)))
assert(isequaln(MKF_SF2.p_rk_g_rkm1, nan(26, 1)))
assert(isequaln(MKF_SF2.p_seq_g_Ykm1, nan(26, 1)))
assert(isequaln(MKF_SF2.xk_est, nan(2, 1)))
assert(isequaln(MKF_SF2.yk_est, nan))

%TODO: Check optional definition with initial states


%% Test convergence to steady-state

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Check steady-state at x0 = [0; 0]
obs = MKF_SF95;
assert(all(abs(obs.filters.Xkp1_est(:, :) - [0; 0]) == 0, [1 2]))
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

% Check steady-state at x0 = [0; 0]
obs = MKF_SF1;
assert(all(abs(obs.filters.Xkp1_est(:, :) - [0; 0]) == 0, [1 2]))
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

% Check steady-state with x0 = [1; 0]
x0 = [1; 0];
io.u_known = u_known;
io.y_meas = true(ny, 1);
nf = 3;
m = 1;
d = 1;
obs = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nf,m, ...
    d,label,x0);
assert(all(abs(obs.filters.Xkp1_est(:, :) - x0) == 0, [1 2]))
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


%% Test sequence updates of MKF_SF_RODD95 on SISO linear system

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Match example 2 from thesis report
x0 = [0; 0];
io.u_known = u_known;
io.y_meas = true(ny, 1);
f = 10;
m = 2;
d = 5;
obs = MKFObserverSF_RODD95(model,io,P0,epsilon,sigma_wp,Q0,R,f,m, ...
    d,label,x0);
assert(strcmp(obs.type, "MKF_SF_RODD95"))
assert(all(abs(obs.filters.Xkp1_est(:, :) - x0) == 0, [1 2]))

% % Generate test simulation data
% nT = 10;
% U_m = zeros(nT+1, sum(u_known));
% % Add a random shock
% Wp = zeros(nT+1, sum(~u_known));
% Wp(5, :) = 1;
% % Compute outputs (n0 measurement noise)
% t = Ts*(0:nT);
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

% Switching sequence
seq = [ ...
    1 1 1 1 1 1 1 1 1 1
    2 1 1 1 1 1 1 1 1 1
    1 1 1 1 1 2 1 1 1 1
    2 1 1 1 1 2 1 1 1 1
];
assert(isequal(cell2mat(obs.seq), seq))

% Initialization of counters
assert(isequal([obs.i obs.i_next], [10 1]))

% Check initial states
assert(isequal(obs.x0, zeros(2, 1)))
assert(isequal(obs.r0, ones(8, 1)))

% Check initialization of probabilities
assert(isequal(obs.rkm1, [0 0 0 0 0 0 0 0]'))
assert(isequal(obs.rk, [1 1 1 1 1 1 1 1]'))
assert(isequaln(obs.p_rk_g_rkm1, nan(8, 1)))
assert(isequaln(obs.p_yk_g_seq_Ykm1, nan(8, 1)))
assert(isequaln(obs.p_seq_g_Ykm1, nan(8, 1)))
assert(isequal(obs.p_seq_g_Yk, ones(8, 1)./8))

% Check initialization of filters
assert(isequal(obs.filters.Pkp1(:,:,1), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,2), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,3), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,4), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,5), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,6), obs.P0))

% Check initial estimates
assert(isequaln(obs.xk_est, nan(2, 1)))
assert(isequaln(obs.yk_est, nan))

% Update at k = 0
i = 1;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [1 2]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 1 2 3 4 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 2 3 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1 1 1 1 1]'))  % obs.r0
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.042 0.042 0.042 0.042 0.042 0.042 0.042 0.042]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.1189 0.0061 0.1189 0.0061 0.1189 0.0061 0.1189 0.0061]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 2 1 2]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4755 0.0245 0.4755 0.0245]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.4755 0.0245 0.4755 0.0245]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 1
i = 2;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [2 3]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 1 1 1]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 2 3 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 2 1 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), [0.9510 0.9510 0.9510 0.9510]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), [0.042 0.042 0.042 0.042]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), [0.4522 0.0233 0.4522 0.0233]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 1 1]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4755 0.0245 0.4755 0.0245]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.4755 0.0245 0.4755 0.0245]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 2
i = 3;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [3 4]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 1 1 1]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 2 3 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1]'))  % obs.r0
assert(isequal(round(obs.p_rk_g_rkm1, 4), [0.9510 0.9510 0.9510 0.9510]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [1.8682 1.0834 1.8682 1.0834]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.4522 0.0233 0.4522 0.0233]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 1 1]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4855 0.0145 0.4855 0.0145]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.4855 0.0145 0.4855 0.0145]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 3
i = 4;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [4 5]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 1 1 1]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 2 3 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), [0.9510 0.9510 0.9510 0.9510]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [2.4676 2.0186 2.4676 2.0186]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.4617 0.0138 0.4617 0.0138]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 1 1]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4881 0.0119 0.4881 0.0119]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.4881 0.0119 0.4881 0.0119]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 4
i = 5;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [5 6]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 1 1 1]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), [0.9510 0.9510 0.9510 0.9510]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [2.8078 2.5333 2.8078 2.5333]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.4641 0.0113 0.4641 0.0113]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 1 1]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4892 0.0108 0.4892 0.0108]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1 1 1 1 1]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.2446 0.2446 0.0054 0.0054 0.2446 0.2446 0.0054 0.0054]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 5
i = 6;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [6 7]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 3 2 4 1 3 2 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 2 3 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1 1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [3.0218 3.0218 2.8436 2.8436 3.0218 3.0218 2.8436 2.8436]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.2326 0.0120 0.0051 0.0003 0.2326 0.0120 0.0051 0.0003]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 2 2]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.9317 0.0193 0.0480 0.0010]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 2 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9317 0.0193 0.0480 0.0010]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 6  *** First non-zero measurement ***
i = 7;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [7 8]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 1 1 1]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 2 3 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 2 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), [0.9510 0.9510 0.9510 0.9510]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.1865 0.2218 0.1865 0.2218]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.8860 0.0184 0.0457 0.0009]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 1 1]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.9281 0.0229 0.0478 0.0012]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.9281 0.0229 0.0478 0.0012]'))

% Check estimates
assert(isequaln(round(obs.xk_est, 4), [0.3719 0.1172]'))
assert(isequaln(round(obs.yk_est, 4), 0.1116))
assert(isequal(round(obs.yk_est - yk, 4), -0.1884))  % -0.1809 MKF_SP_RODD

% Update at k = 7
i = 8;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [8 9]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 1 1 1]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 2 3 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), [0.9510 0.9510 0.9510 0.9510]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0166 0.0335 0.5807 0.6228]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.8826 0.0218 0.0455 0.0011]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 1 1]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.3450 0.0172 0.6214 0.0165]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.3450 0.0172 0.6214 0.0165]'))

% Check estimates
assert(isequaln(round(obs.xk_est, 4), [1.3001 0.8975]'))
assert(isequaln(round(obs.yk_est, 4), 0.3900))
assert(isequal(round(obs.yk_est - yk, 4), -0.1200))  % -0.0563 MKF_SP_RODD

% Update at k = 8
i = 9;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [9 10]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 1 1 1]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 2 3 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), [0.9510 0.9510 0.9510 0.9510]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0084 0.0242 1.9563 1.9729]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.3281 0.0163 0.5909 0.0156]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 1 1]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.0023 0.0003 0.9714 0.0259]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.0023 0.0003 0.9714 0.0259]'))

% Check estimates
assert(isequaln(round(obs.xk_est, 4), [2.2343 1.1866]'))
assert(isequaln(round(obs.yk_est, 4), 0.6703))
assert(isequal(round(obs.yk_est - yk, 4), 0.0133))  % -0.0101 MKF_SP_RODD

% Update at k = 9
i = 10;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [10 1]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 1 1 1]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), [0.9510 0.9510 0.9510 0.9510]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0114 0.0369 2.3413 2.3806]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.0022 0.0003 0.9238 0.0247]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 1 1]'))
assert(isequal(obs.merged.rk, seq(:, i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.0000 0.0000 0.9736 0.0264]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1 1 1 1 1]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.0000 0.0000 0.0000 0.0000 0.4868 0.4868 0.0132 0.0132]'))

% Check estimates
assert(isequaln(round(obs.xk_est, 4), [2.6247 1.1289]'))
assert(isequaln(round(obs.yk_est, 4), 0.7874))
assert(isequal(round(obs.yk_est - yk, 4), 0.0275))


%% Test sequence updates of MKF_SF_RODD on SISO linear system

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Match example 2 from thesis report
x0 = [0; 0];
io.u_known = u_known;
io.y_meas = true(ny, 1);
nf = 2;
m = 2;
d = 5;
obs = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nf,m, ...
    d,"MKF_SF_RODD",x0);
assert(strcmp(obs.type, "MKF_SF_RODD"))
assert(all(abs(obs.filters.Xkp1_est(:, :) - x0) == 0, [1 2]))

% % Generate test simulation data
% nT = 10;
% U_m = zeros(nT+1, sum(u_known));
% % Add a random shock
% Wp = zeros(nT+1, sum(~u_known));
% Wp(5, :) = 1;
% % Compute outputs (n0 measurement noise)
% t = Ts*(0:nT);
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

% Switching sequence
seq = [ ...
    1 1
    2 1
    1 2
    2 2
];
assert(isequal(cell2mat(obs.seq), seq))

% Initialization of counters
assert(isequal([obs.i obs.i_next], [2 1]))

% Check initial states
assert(isequal(obs.x0, zeros(2, 1)))
assert(isequal(obs.r0, [1 1 1 1 2 2 2 2]'))

% Check initialization of probabilities
assert(isequal(obs.rkm1, [0 0 0 0 0 0 0 0]'))
assert(isequal(obs.rk, [1 1 1 1 2 2 2 2]'))
assert(isequaln(obs.p_rk_g_rkm1, nan(8, 1)))
assert(isequaln(obs.p_yk_g_seq_Ykm1, nan(8, 1)))
assert(isequaln(obs.p_seq_g_Ykm1, nan(8, 1)))
assert(isequal(obs.p_seq_g_Yk, ones(8, 1)./8))

% Check initialization of filters
assert(isequal(obs.filters.Pkp1(:,:,1), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,2), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,3), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,4), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,5), obs.P0))
assert(isequal(obs.filters.Pkp1(:,:,6), obs.P0))

% Check initial estimates
assert(isequaln(obs.xk_est, nan(2, 1)))
assert(isequaln(obs.yk_est, nan))

% Update at k = 0  % no branching and merging
i = 1;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [1 2]))
assert(isequal([obs.id obs.id_next], [1 2]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 1 2 3 4 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1 2 2 2 2]'))  % obs.r0
assert(isequaln(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
%assert(isequal( ...
%    obs.p_rk_g_rkm1, ...
%    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
%))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.042 0.042 0.042 0.042 0.042 0.042 0.042 0.042]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.1189 0.0061 0.1189 0.0061 0.1189 0.0061 0.1189 0.0061]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 2 1 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4755 0.0245 0.4755 0.0245]'))

%TODO: Need to check above and change below after making change 2022-12-07

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.2377 0.0123 0.2377 0.0123 0.2377 0.0123 0.2377 0.0123]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 1
i = 2;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [1 2]))
assert(isequal([obs.id obs.id_next], [2 3]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 1 2 3 4 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1 2 2 2 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.042 0.042 0.042 0.042 0.042 0.042 0.042 0.042]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.2377 0.0123 0.2377 0.0123 0.2377 0.0123 0.2377 0.0123]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 2 1 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4755 0.0245 0.4755 0.0245]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.2377 0.0123 0.2377 0.0123 0.2377 0.0123 0.2377 0.0123]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 2
i = 3;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [1 2]))
assert(isequal([obs.id obs.id_next], [3 4]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 1 2 3 4 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1 2 2 2 2]'))  % obs.r0
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [1.8682 1.5820 1.8682 1.5820 1.8682 1.5820 1.8682 1.5820]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.2377 0.0123 0.2377 0.0123 0.2377 0.0123 0.2377 0.0123]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 2 1 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4791 0.0209 0.4791 0.0209]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.2395 0.0105 0.2395 0.0105 0.2395 0.0105 0.2395 0.0105]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 3
i = 4;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [1 2]))
assert(isequal([obs.id obs.id_next], [4 5]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 1 2 3 4 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 1 1 2 2 2 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [2.4676 1.7906 2.4676 1.7906 2.4676 1.7906 2.4676 1.7906]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.2395 0.0105 0.2395 0.0105 0.2395 0.0105 0.2395 0.0105]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 2 1 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4847 0.0153 0.4847 0.0153]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.2423 0.0077 0.2423 0.0077 0.2423 0.0077 0.2423 0.0077]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 4
i = 5;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [1 2]))
assert(isequal([obs.id obs.id_next], [5 1]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 2 1 2 3 4 3 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [2.8078 1.8085 2.8078 1.8085 2.8078 1.8085 2.8078 1.8085]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.2423 0.0077 0.2423 0.0077 0.2423 0.0077 0.2423 0.0077]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 2 1 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.4900 0.0100 0.4900 0.0100]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 2 2 1 1 2 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.2450 0.2450 0.0050 0.0050 0.2450 0.2450 0.0050 0.0050]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 5  % Into next detection interval
i = 6;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [2 1]))
assert(isequal([obs.id obs.id_next], [1 2]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 3 2 4 1 3 2 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 2 2 1 1 2 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [3.0218 3.0218 1.8086 1.8086 3.0218 3.0218 1.8086 1.8086]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.2330 0.0120 0.0048 0.0002 0.2330 0.0120 0.0048 0.0002]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 2 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.9395 0.0115 0.0484 0.0006]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.4698 0.0242 0.0057 0.0003 0.4698 0.0242 0.0057 0.0003]'))

% Check estimates
assert(isequaln(obs.xk_est, [0; 0]))
assert(isequaln(obs.yk_est, 0))
assert(isequal(obs.yk_est - yk, 0))

% Update at k = 6  *** First non-zero measurement ***
i = 7;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [2 1]))
assert(isequal([obs.id obs.id_next], [2 3]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 3 2 4 1 3 2 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 2 2 1 1 2 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.1865 0.1865 0.7172 0.7172 0.1865 0.1865 0.7172 0.7172]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.4698 0.0242 0.0057 0.0003 0.4698 0.0242 0.0057 0.0003]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 2 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.9083 0.0427 0.0468 0.0022]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.4542 0.0234 0.0213 0.0011 0.4542 0.0234 0.0213 0.0011]'))

% Check estimates
assert(isequaln(round(obs.xk_est, 4), [0.3898 0.1385]'))
assert(isequaln(round(obs.yk_est, 4), 0.1169))
assert(isequal(round(obs.yk_est - yk, 4), -0.1831))  % -0.1809 MKF_SP_RODD

% Update at k = 7
i = 8;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [2 1]))
assert(isequal([obs.id obs.id_next], [3 4]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 3 2 4 1 3 2 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 2 2 1 1 2 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0166 0.2003 1.4956 1.3873 0.0166 0.2003 1.4956 1.3873]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.4542 0.0234 0.0213 0.0011 0.4542 0.0234 0.0213 0.0011]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 2 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.1653 0.6986 0.1027 0.0334]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.0826 0.0514 0.3493 0.0167 0.0826 0.0514 0.3493 0.0167]'))

% Check estimates
assert(isequaln(round(obs.xk_est, 4), [1.3865 0.7297]'))
assert(isequaln(round(obs.yk_est, 4), 0.4159))
assert(isequal(round(obs.yk_est - yk, 4), -0.0941))  % -0.0563 MKF_SP_RODD

% Update at k = 8
i = 9;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [2 1]))
assert(isequal([obs.id obs.id_next], [4 5]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 3 2 4 1 3 2 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 1 2 2 1 1 2 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0084 1.5880 2.1922 1.7750 0.0084 1.5880 2.1922 1.7750]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.0826 0.0514 0.3493 0.0167 0.0826 0.0514 0.3493 0.0167]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 2 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.0008 0.8725 0.0929 0.0338]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.0004 0.0465 0.4363 0.0169 0.0004 0.0465 0.4363 0.0169]'))

% Check estimates
assert(isequaln(round(obs.xk_est, 4), [2.0598 0.9155]'))
assert(isequaln(round(obs.yk_est, 4), 0.6179))
assert(isequal(round(obs.yk_est - yk, 4), -0.0391))  % -0.0101 MKF_SP_RODD

% Update at k = 9
i = 10;
uk = U_m(i,:)';
yk = Y_m(i,:)';
obs.update(yk, uk);
assert(isequal([obs.i obs.i_next], [2 1]))
assert(isequal([obs.id obs.id_next], [5 1]))

% Check merging and branching indices
assert(isequal(obs.idx_modes{obs.i}, [1 2 1 2 1 2 1 2]'))
assert(isequal(obs.idx_merge{obs.i}, [1 3 2 4 1 3 2 4]'))
assert(isequal(obs.idx_branch{obs.i_next}, [1 1 2 2 3 3 4 4]'))

% Check branched hypotheses probabilities
assert(isequal(obs.rkm1, [1 2 1 2 1 2 1 2]'))
assert(isequal(round(obs.p_rk_g_rkm1, 4), ...
    [0.9510 0.0490 0.9510 0.0490 0.9510 0.0490 0.9510 0.0490]'))
assert(isequal( ...
    obs.p_rk_g_rkm1, ...
    prob_rk(obs.idx_modes{obs.i}, obs.T(obs.rkm1, :)') ...
))
assert(isequal(round(obs.p_yk_g_seq_Ykm1, 4), ...
    [0.0114 1.8008 2.6291 1.8089 0.0114 1.8008 2.6291 1.8089]'))
assert(isequal(round(obs.p_seq_g_Ykm1, 4), ...
    [0.0004 0.0465 0.4363 0.0169 0.0004 0.0465 0.4363 0.0169]'))

% Check merged hypotheses
assert(isequal(obs.merged.rk, [1 1 2 2]'))
assert(isequal(obs.merged.rk, seq(:, obs.i)))
assert(isequal(round(obs.merged.p_seq_g_Yk, 4), ...
    [0.0000 0.9094 0.0663 0.0242]'))

% Branched hypotheses at end of this time instant
assert(isequal(obs.rk, [1 1 1 1 2 2 2 2]'))
assert(isequal(round(obs.p_seq_g_Yk, 4), ...
    [0.0000 0.0000 0.4547 0.4547 0.0332 0.0332 0.0121 0.0121]'))

% Check estimates
assert(isequaln(round(obs.xk_est, 4), [2.4327 0.9387]'))
assert(isequaln(round(obs.yk_est, 4),  0.7298))
assert(isequal(round(obs.yk_est - yk, 4), -0.0301))


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
observers = {KF2, KF3, SKF, MKF3, MKF4, MKF_SF95, MKF_SF1, MKF_SF2};
% Note: KF1 is too slow to pass static error test here
% Results of last observer simulation are plotted.

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

    % Compute mean-squared error
    Y_est = sim_results.Y_est;
    MSE.(obs.label) = mean((Y_est - Y).^2);

    % Save updated observer
    observers{i} = obs;

end

MSE_test_values = struct(...
    'MKF_SF95', 0.000764, ...
    'MKF_SF1', 0.000448, ...  % c.f. 0.000491 MKF_SP1
    'MKF_SF2', 0.000565, ...  % c.f. 0.000492 MKF_SP2
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

% X_est = sim_results.X_est;
% E_obs = sim_results.E_obs;
% K_obs = sim_results.K_obs;
% trP_obs = sim_results.trP_obs;

% % Display results of last simulation
% table(t,alpha,U,P,Wp,X,Y,Y_m,X_est,Y_est,E_obs)

% % Display gains and trace of covariance matrix
% table(t, cell2mat(K_obs), trP_obs, ...
%     'VariableNames',{'t', 'K{1}, K{2}', 'trace(P)'})

% % Plot of inputs and outputs
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

% % Plot of merged hypothesis probabilities
% switch obs.type
%     case {'MKF_SF', 'MKF_SF_RODD'}
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

% % Plot of trace of filter covariance matrices
% switch obs.type
%     case {'MKF', 'MKF_S', 'MKF_SP', 'MKF_SP_RODD', ...
%             'MKF_SF', 'MKF_SF_RODD', 'MKF_SF_RODD95'}
%         figure(12); clf
%         t = Ts*(0:nT)';
%         ax_labels = {'$t$', 'MKF filter ($\Gamma(k)$)', '$Tr(P(k))$'};
%         make_waterfall_plot(t, sim_results.trP_obs_f, [0 1], ax_labels, [0 82]);
%         filename = sprintf('rod_mkf_observer_test_trP_wfplot.pdf');
%         save_fig_to_pdf(fullfile(plot_dir, filename));
%         title('Trace of covariance matrices')
% end


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

% Multiple model observer with sequence fusion 1
label = "MKF_SF1";
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
nf = 3;  % Number of detection intervals in fusion horizon
m = 1;  % Maximum number of shocks over fusion horizon
d = 5;  % Length of detection intervals in sample periods
io.u_known = u_known;
io.y_meas = true(ny, 1);
% TODO: Replace cell2mat(sigma_wp') with sigma_wp when MKFObserverSF_RODD
% updated
MKF_SF1 = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nf,m,d,label);

% Multiple model observer with sequence fusion 2
label = "MKF_SF2";
P0 = 1000*eye(n);
Q0 = diag([0.01 0.01 0 0]);
R = diag(sigma_M.^2);
nf = 5;  % Number of detection intervals in fusion horizon
m = 2;  % Maximum number of shocks over fusion horizon
d = 3;  % Length of detection intervals in sample periods
io.u_known = u_known;
io.y_meas = true(ny, 1);
MKF_SF2 = MKFObserverSF_RODD(model,io,P0,epsilon,sigma_wp,Q0,R,nf,m,d,label);

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
observers = {KF3, MKF_SF1, MKF_SF2, MKF3, MKF4, SKF};

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

% Results on 2022-11-29 after modifying MKF_SF_DI
MSE_test_values = struct(...
    'KF3', [0.000296 0.000433], ...
    'MKF_SF1', [0.000870 0.000747], ...  % c.f. [0.000308 0.000758] MKF_SP1
    'MKF_SF2', [0.000265 0.000305], ...  % c.f. [0.000322 0.000773] MKF_SP2
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
