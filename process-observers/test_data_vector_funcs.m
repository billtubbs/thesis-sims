% Test the following functions used to run observers
% in Simulink
% 
%  - make_data_vectors.m
%  - unpack_data_vectors.m
%  - get_obs_vars.m
%  - get_obs_vars_vecs.m
%  - set_obs_vars_vecs.m
%


%% Test make_data_vectors

vdata = make_data_vectors({0});
assert(isequal(vdata.vecs, {0}))
assert(isequal(vdata.types, {'double'}))
assert(isequal(vdata.dims, {[1 1]}))
assert(isequal(vdata.n_els, {[1]}))

vdata = make_data_vectors({1, 2, 3});
assert(isequal(vdata.vecs, {1, 2, 3}))
assert(isequal(vdata.types, {'double', 'double', 'double'}))
assert(isequal(vdata.dims, {[1 1], [1 1], [1 1]}))
assert(isequal(vdata.n_els, {1, 1, 1}))

vdata = make_data_vectors({1, [2; 3], [4 6 8; 5 7 9]});
assert(isequal(vdata.vecs, {1, [2 3], [4 5 6 7 8 9]}))
assert(isequal(vdata.types, {'double', 'double', 'double'}))
assert(isequal(vdata.dims, {[1 1], [2 1], [2 3]}))
assert(isequal(vdata.n_els, {1, 2, 6}))

x1 = {1 3 5; 2 4 6};
vdata_x1 = make_data_vectors({x1});
assert(isequal(vdata_x1.vecs, {1:6}))
assert(isequal(vdata_x1.types, {{'double', 'double', 'double'; 'double', 'double', 'double'}}))
assert(isequal(vdata_x1.dims, {{[1 1], [1 1], [1 1]; [1 1], [1 1], [1 1]}}))
assert(isequal(vdata_x1.n_els, {6}))

x2 = {[7 8 9], 10, [11; 12], [13 15 17; 14 16 18]};
vdata_x2 = make_data_vectors({x2});
assert(isequal(vdata_x2.vecs, {7:18}))
assert(isequal(vdata_x2.types, {{'double', 'double', 'double', 'double'}}))
assert(isequal(vdata_x2.dims, {{[1 3], [1 1], [2 1], [2 3]}}))
assert(isequal(vdata_x2.n_els, {12}))

vdata_y = make_data_vectors({x1, x2, [19; 20; 21]});
assert(isequal(cell2mat(vdata_y.vecs), 1:21))
assert(isequal(vdata_y.types, {vdata_x1.types{1}, vdata_x2.types{1}, 'double'}))
assert(isequal(vdata_y.dims, {vdata_x1.dims{1}, vdata_x2.dims{1}, [3 1]}))
assert(isequal(vdata_y.n_els, {6, 12, 3}))

A = [1 3 5; 2 4 6];
vdata = make_data_vectors({A, x2, 19, [20; 21]});
assert(isequal(cell2mat(vdata.vecs), 1:21))
assert(isequal(vdata.types, {'double', vdata_x2.types{1}, 'double', 'double'}))
assert(isequal(vdata.dims, {size(A), vdata_x2.dims{1}, [1 1], [2 1]}))
assert(isequal(vdata.n_els, {6, 12, 1, 2}))


%% Test unpack_data_vectors

vdata = make_data_vectors({1});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {1}))

vdata = make_data_vectors({1, 2});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {1, 2}))

vdata = make_data_vectors({1, [2; 3], [4 5]});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {1, [2; 3], [4 5]}))

% Cell array
a = 1;
b = [2 3; 4 5];
c = {6, [7; 8], 9};
vdata = make_data_vectors({a, b, c});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {a, b, c}))

% Cell array
a = 1;
b = [2 3; 4 5];
c = {6, [7; 8], 9};
d = {10; 11};
vdata = make_data_vectors({a b; c d});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {a b; c d}))

% Nested cell arrays
a = 1;
b = [2 3; 4 5];
c = {6, {[7 8 9], [10; 11]}, 12};
vdata = make_data_vectors({a, b, c});
vars = unpack_data_vectors(vdata);
assert(isequal(vars, {a, b, c}))


%% Test make_data_vectors with different numeric data types

i = int16(10);
vdata = make_data_vectors({i}, 'int16');
assert(isequal(vdata.vecs, {int16(10)}))
assert(isequal(vdata.types, {'int16'}))
assert(isequal(vdata.dims, {[1 1]}))

vdata = make_data_vectors({uint8(1), uint8([2; 3]), ...
    uint8([4 6 8; 5 7 9])}, 'uint8');
assert(isequal(vdata.vecs, {uint8(1), uint8([2 3]), uint8([4 5 6 7 8 9])}))
assert(isequal(vdata.types, {'uint8', 'uint8', 'uint8'}))
assert(isequal(vdata.dims, {[1 1], [2 1], [2 3]}))


%% Test get_obs_vars.m, get_obs_vars_vecs.m and set_obs_vars_vecs.m

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

vars = get_obs_vars(KFPSS1);
assert(isequal(fieldnames(vars), {'xkp1_est', 'ykp1_est'}'))
assert(isequal(vars.xkp1_est, KFPSS1.xkp1_est))
assert(isequal(vars.ykp1_est, KFPSS1.ykp1_est))
vars_vecs = get_obs_vars_vecs(KFPSS1);
assert(isequal(vars_vecs, [KFPSS1.xkp1_est' KFPSS1.ykp1_est']))
obs_copy = KFPSS1.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vars_vecs');
% TODO: Why is vars_vecs a column vector here?
% TODO: Need to test it actually set the variables correctly

vars = get_obs_vars(KFFSS1);
assert(isequal(fieldnames(vars), {'xkp1_est'}'))
assert(isequal(vars.xkp1_est, KFFSS1.xkp1_est))
vars_vecs = get_obs_vars_vecs(KFFSS1);
assert(isequal(vars_vecs, KFFSS1.xkp1_est'))
obs_copy = KFFSS1.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vars_vecs');

vars = get_obs_vars(LB1);
assert(isequal(fieldnames(vars), {'xkp1_est', 'ykp1_est'}'))
assert(isequal(vars.xkp1_est, LB1.xkp1_est))
assert(isequal(vars.ykp1_est, LB1.ykp1_est))
vars_vecs = get_obs_vars_vecs(LB1);
assert(isequal(vars_vecs, [LB1.xkp1_est' LB1.ykp1_est']))
obs_copy = LB1.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vars_vecs');

vars = get_obs_vars(KF1);
assert(isequal(fieldnames(vars), {'xkp1_est', 'Pkp1'}'))
assert(isequal(vars.xkp1_est, KF1.xkp1_est))
assert(isequal(vars.Pkp1, KF1.Pkp1))
vars_vecs = get_obs_vars_vecs(KF1);
assert(isequal(vars_vecs, [KF1.xkp1_est' KF1.Pkp1(:)']))
obs_copy = KF1.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vars_vecs');

vars = get_obs_vars(KFP3);
assert(isequal(fieldnames(vars), {'xkp1_est', 'ykp1_est', 'Pkp1'}'))
assert(isequal(vars.xkp1_est, KFP3.xkp1_est))
assert(isequal(vars.ykp1_est, KFP3.ykp1_est))
assert(isequal(vars.Pkp1, KFP3.Pkp1))
vars_vecs = get_obs_vars_vecs(KFP3);
assert(isequal(vars_vecs, [KFP3.xkp1_est' KFP3.ykp1_est' KFP3.Pkp1(:)']))
obs_copy = KFP3.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vars_vecs');

% vars = get_obs_vars(MKF_SF1);
% assert(isequal(fieldnames(vars), {'xkp1_est', 'p_seq_g_Yk', ...
%     'xkp1_est_f', 'Pkp1_f', 'int16'}'))
% assert(isequal(vars.xkp1_est, MKF_SF1.xkp1_est))
% assert(isequal(vars.p_seq_g_Yk, MKF_SF1.p_seq_g_Yk))
% assert(isequal(vars.xkp1_est_f, MKF_SF1.filters.Xkp1_est))
% assert(isequal(vars.Pkp1_f, MKF_SF1.filters.Pkp1))
% assert(isequal(fieldnames(vars.int16), {'rk', 'i', 'i_next', 'i2', 'i2_next'}'))
% assert(isequal(vars.int16.rk, MKF_SF1.rk))
% assert(isequal(vars.int16.i, MKF_SF1.i))
% assert(isequal(vars.int16.i_next, MKF_SF1.i_next))
% assert(isequal(vars.int16.i2, MKF_SF1.i2))
% assert(isequal(vars.int16.i2_next, MKF_SF1.i2_next))
% [vec_double, vec_int16] = get_obs_vars_vecs(MKF_SF1);
% assert(isequal(vec_double, [ ...
%     MKF_SF1.xkp1_est' MKF_SF1.p_seq_g_Yk' ...
%     reshape(MKF_SF1.filters.Xkp1_est, 1, MKF_SF1.nh*MKF_SF1.n) ...
%     reshape(MKF_SF1.filters.Pkp1,1,MKF_SF1.nh*MKF_SF1.n^2) ...
% ]))
% assert(isequal(vec_int16, [ ...
%     MKF_SF1.rk' MKF_SF1.i MKF_SF1.i_next ...
%     MKF_SF1.i2 MKF_SF1.i2_next ...
% ]))
% obs_copy = MKF_SF1.copy();
% obs_copy = set_obs_vars_vecs(obs_copy, vec_double', vec_int16');

vars = get_obs_vars(MKF_SP1);
assert(isequal(fieldnames(vars), {'xkp1_est', 'p_seq_g_Yk', ...
    'xkp1_est_f', 'Pkp1_f', 'int16'}'))
assert(isequal(vars.xkp1_est, MKF_SP1.xkp1_est))
assert(isequal(vars.p_seq_g_Yk, MKF_SP1.p_seq_g_Yk))
assert(isequal(vars.xkp1_est_f, MKF_SP1.filters.Xkp1_est))
assert(isequal(vars.Pkp1_f, MKF_SP1.filters.Pkp1))
assert(isequal(fieldnames(vars.int16), {'rk', 'id', 'id_next', ...
    'f_main', 'f_hold'}'))
assert(isequal(vars.int16.rk, MKF_SP1.rk))
assert(isequal(vars.int16.id, MKF_SP1.id))
assert(isequal(vars.int16.id_next, MKF_SP1.id_next))
assert(isequal(vars.int16.f_main, MKF_SP1.f_main))
assert(isequal(vars.int16.f_hold, MKF_SP1.f_hold))
[vec_double, vec_int16] = get_obs_vars_vecs(MKF_SP1);
assert(isequal(vec_double, [ ...
    MKF_SP1.xkp1_est' MKF_SP1.p_seq_g_Yk' ...
    reshape(MKF_SP1.filters.Xkp1_est, 1, MKF_SP1.nh*MKF_SP1.n) ...
    reshape(MKF_SP1.filters.Pkp1,1,MKF_SP1.nh*MKF_SP1.n^2) ...
]))
assert(isequal(vec_int16, [MKF_SP1.rk' MKF_SP1.id ...
    MKF_SP1.id_next MKF_SP1.f_main MKF_SP1.f_hold]))
obs_copy = MKF_SP1.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vec_double', vec_int16');

% Test switching Kalman filter
SKF = SKFObserver(obs_models,P0,"SKF");
vars = get_obs_vars(SKF);
assert(isequal(fieldnames(vars), {'xkp1_est', 'Pkp1'}'))
assert(isequal(vars.xkp1_est, KF1.xkp1_est))
assert(isequal(vars.Pkp1, KF1.Pkp1))
vars_vecs = get_obs_vars_vecs(SKF);
assert(isequal(vars_vecs, [SKF.xkp1_est' SKF.Pkp1(:)']))
obs_copy = SKF.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vars_vecs');

% Test switching Kalman filter - scheduled
seq = ones(1, 11);
SKF_S = SKFObserverS(obs_models,P0,seq,"SKF_S");
vars = get_obs_vars(SKF_S);
assert(isequal(fieldnames(vars), {'xkp1_est', 'Pkp1', 'int16'}'))
assert(isequal(vars.xkp1_est, KF1.xkp1_est))
assert(isequal(vars.Pkp1, KF1.Pkp1))
assert(isequal(fieldnames(vars.int16), {'i', 'i_next'}'))
assert(isequal(vars.int16.i, SKF_S.i))
assert(isequal(vars.int16.i_next, SKF_S.i_next))
[vec_double, vec_int16] = get_obs_vars_vecs(SKF_S);
assert(isequal(vec_double, [SKF_S.xkp1_est' SKF_S.Pkp1(:)']))
assert(isequal(vec_int16, [SKF_S.i SKF_S.i_next]))
obs_copy = SKF_S.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vec_double', vec_int16');

% Test basic multi-model observer classes
T = [0.95 0.05; 0.2 0.8];
P0 = [100 1000] .* eye(2);
r0 = [1 1 2]';
MKF1 = MKFObserver(obs_models,P0,T,r0,"MKF1");
vars = get_obs_vars(MKF1);
assert(isequal(fieldnames(vars), ...
    {'xkp1_est', 'p_seq_g_Yk', 'xkp1_est_f', 'Pkp1_f', 'int16'}' ...
))
assert(isequal(vars.xkp1_est, MKF1.xkp1_est))
assert(isequal(vars.p_seq_g_Yk, MKF1.p_seq_g_Yk))
assert(isequal(vars.xkp1_est_f, MKF1.filters.Xkp1_est))
assert(isequal(vars.Pkp1_f, MKF1.filters.Pkp1))
assert(isequal(fieldnames(vars.int16), {'rk'}'))
assert(isequal(vars.int16.rk, MKF1.rk))
[vec_double, vec_int16] = get_obs_vars_vecs(MKF1);
assert(isequal(vec_double, ...
    [MKF1.xkp1_est' MKF1.p_seq_g_Yk' MKF1.filters.Xkp1_est(:)' MKF1.filters.Pkp1(:)'] ...
))
assert(isequal(vec_int16, MKF1.rk'))
obs_copy = MKF1.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vec_double', vec_int16');

% MKF observer with switching schedule
T = [0.95 0.05; 0.2 0.8];
P0 = [100 1000] .* eye(2);
seq = {ones(1, 11); 2*ones(1, 11); randi(2, 1, 11)};
MKF2 = MKFObserverS(obs_models,P0,seq,T,"MKF1");
vars = get_obs_vars(MKF2);
assert(isequal(fieldnames(vars), ...
    {'xkp1_est', 'p_seq_g_Yk', 'xkp1_est_f', 'Pkp1_f', 'int16'}' ...
))
assert(isequal(vars.xkp1_est, MKF2.xkp1_est))
assert(isequal(vars.p_seq_g_Yk, MKF2.p_seq_g_Yk))
assert(isequal(vars.xkp1_est_f, MKF2.filters.Xkp1_est))
assert(isequal(vars.Pkp1_f, MKF2.filters.Pkp1))
assert(isequal(fieldnames(vars.int16), {'i', 'i_next'}'))
assert(isequal(vars.int16.i, MKF2.i))
assert(isequal(vars.int16.i_next, MKF2.i_next))
[vec_double, vec_int16] = get_obs_vars_vecs(MKF2);
assert(isequal(vec_double, ...
    [MKF2.xkp1_est' MKF2.p_seq_g_Yk' MKF2.filters.Xkp1_est(:)' MKF2.filters.Pkp1(:)'] ...
))
assert(isequal(vec_int16, [MKF2.i MKF2.i_next]))
obs_copy = MKF2.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vec_double', vec_int16');

% MKF observer with sequence pruning
nh = 5;
n_min = 3;
MKF_SP = MKFObserverSP(obs_models,P0,T,nh,n_min,"MKF_SP2");
vars = get_obs_vars(MKF_SP);
assert(isequal(fieldnames(vars), {'xkp1_est', 'p_seq_g_Yk', ...
    'xkp1_est_f', 'Pkp1_f', 'int16'}'))
assert(isequal(vars.xkp1_est, MKF_SP.xkp1_est))
assert(isequal(vars.p_seq_g_Yk, MKF_SP.p_seq_g_Yk))
assert(isequal(vars.xkp1_est_f, MKF_SP.filters.Xkp1_est))
assert(isequal(vars.Pkp1_f, MKF_SP.filters.Pkp1))
assert(isequal(fieldnames(vars.int16), {'rk', 'f_main', 'f_hold'}'))
assert(isequal(vars.int16.rk, MKF_SP.rk))
assert(isequal(vars.int16.f_main, MKF_SP.f_main))
assert(isequal(vars.int16.f_hold, MKF_SP.f_hold))
[vec_double, vec_int16] = get_obs_vars_vecs(MKF_SP);
assert(isequal(vec_double, [ ...
    MKF_SP.xkp1_est' MKF_SP.p_seq_g_Yk' ...
    reshape(MKF_SP.filters.Xkp1_est, 1, MKF_SP.nh*MKF_SP1.n) ...
    reshape(MKF_SP.filters.Pkp1, 1, MKF_SP.nh*MKF_SP.n^2) ...
]))
assert(isequal(vec_int16, [MKF_SP.rk' MKF_SP.f_main MKF_SP.f_hold]))
obs_copy = MKF_SP.copy();
obs_copy = set_obs_vars_vecs(obs_copy, vec_double', vec_int16');


%% Test in simulation with Kalman filter object

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_wp{1}(2), sigma_wp{1}(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
Y = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Make copies of KF1
obs1 = KF1.copy();
obs2 = KF1.copy();

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vdata = make_data_vectors(struct2cell(vars)');

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT

    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1.update(yk, uk);

    % Save observer estimates
    Xk_est{1}(i+1,:) = obs1.xk_est';
    Yk_est{1}(i+1,:) = obs1.yk_est';

    % Unpack vdata struct and reconstruct observer from KF1
    vars = unpack_data_vectors(vdata);
    assert(isequal(vars{1}, obs2.xkp1_est))
    assert(isequal(vars{2}, obs2.Pkp1))
    obs2 = KF1;  % makes a new copy
    obs2.xkp1_est = vars{1};
    obs2.Pkp1 = vars{2};

    % Observer updates
    obs2.update(yk, uk);

    % Convert dynamic variables to vdata struct
    vars = get_obs_vars(obs2);
    vdata = make_data_vectors(struct2cell(vars)');

    % Save observer estimates
    Xk_est{2}(i+1,:) = obs2.xk_est';
    Yk_est{2}(i+1,:) = obs2.yk_est';

end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))


%% Test with MKF_SP1 observer

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_wp{1}(2), sigma_wp{1}(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
Y = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Make copies of MKF_SP1
obs1 = MKF_SP1.copy();
obs2 = MKF_SP1.copy();

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vars_double = {vars.xkp1_est, vars.p_seq_g_Yk, ...
    vars.xkp1_est_f, vars.Pkp1_f};
vdata = make_data_vectors(vars_double);
vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT

    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1.update(yk, uk);

    % Save observer estimates
    Xk_est{1}(i,:) = obs1.xk_est';
    Yk_est{1}(i,:) = obs1.yk_est';

    % Unpack vdata struct and copy values back to observer
    vars_double = unpack_data_vectors(vdata);
    vars_int16 = unpack_data_vectors(vdata_int16);
    vars = struct();
    vars.xkp1_est = vars_double{1};
    %vars.Pkp1 = vars_double{2};  % do we need this?
    vars.p_seq_g_Yk = vars_double{2};
    vars.xkp1_est_f = vars_double{3};
    vars.Pkp1_f = vars_double{4};
    vars.int16.rk = vars_int16{1};
    vars.int16.id = vars_int16{2};
    vars.int16.id_next = vars_int16{3};
    vars.int16.f_main = vars_int16{4};
    vars.int16.f_hold = vars_int16{5};

    % Build a copy from the variable data
    obs2_restored = set_obs_vars(MKF_SP1.copy(), vars);
    assert(isequal(obs2_restored.xkp1_est, obs2.xkp1_est))
    %assert(isequal(obs2_restored.Pkp1, obs2.Pkp1))
    assert(isequal(obs2_restored.p_seq_g_Yk, obs2.p_seq_g_Yk))
    assert(isequal(obs2_restored.filters.Xkp1_est, obs2.filters.Xkp1_est))
    assert(isequal(obs2_restored.filters.Pkp1, obs2.filters.Pkp1))
    assert(isequal(obs2_restored.rk, obs2.rk))
    assert(isequal(obs2_restored.f_main, obs2.f_main))
    assert(isequal(obs2_restored.f_hold, obs2.f_hold))

    % Observer updates
    obs2.update(yk, uk);

    % Re-convert dynamic variables to vdata struct for next
    % time step
    vars = get_obs_vars(obs2);
    vars_double = {vars.xkp1_est, vars.p_seq_g_Yk, ...
        vars.xkp1_est_f, vars.Pkp1_f};
    vdata = make_data_vectors(vars_double);
    vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');

    % Save observer estimates
    Xk_est{2}(i,:) = obs2.xk_est';
    Yk_est{2}(i,:) = obs2.yk_est';

end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))


%% Test with MKF_SF1 observer

% Load system and disturbance model from file
sys_rodin_step

% Load observers from file
obs_rodin_step

% Sequenece length
nT = 20;

% Simulate system
X0 = zeros(n,1);
t = Ts*(0:nT)';
Wp = sample_random_shocks(nT+1, epsilon, sigma_wp{1}(2), sigma_wp{1}(1));
U = zeros(nT+1,1);
U(t>=5) = 1;
[Y,T,X] = lsim(Gpss,[U Wp],t,X0);
V = sigma_M*randn(nT+1, 1);
Ym = Y + V;  % measurement

% Make copies of MKF_SF1
obs1 = MKF_SF1.copy();
obs2 = MKF_SF1.copy();

% Convert dynamic variables to vdata struct
vars = get_obs_vars(obs2);
vars_double = {vars.xkp1_est, vars.p_seq_g_Yk, ...
    vars.xkp1_est_f, vars.Pkp1_f};
vdata = make_data_vectors(vars_double);
vdata_int16 = make_data_vectors(struct2cell(vars.int16), 'int16');

% Simulate observer
Xk_est = {nan(nT+1,n), nan(nT+1,n)};
Yk_est = {nan(nT+1,ny), nan(nT+1,ny)};
for i = 1:nT

    uk = U(i,:)';
    yk = Ym(i,:)';
    obs1.update(yk, uk);

    % Save observer estimates
    Xk_est{1}(i+1,:) = obs1.xk_est';
    Yk_est{1}(i+1,:) = obs1.yk_est';

    % Unpack vdata struct and copy values back to observer
    vars_double = unpack_data_vectors(vdata);
    vars_int16 = unpack_data_vectors(vdata_int16);
    vars = struct();
    vars.xkp1_est = vars_double{1};
    %vars.Pkp1 = vars_double{2};
    vars.p_seq_g_Yk = vars_double{2};
    vars.xkp1_est_f = vars_double{3};
    vars.Pkp1_f = vars_double{4};
    vars.int16.rk = vars_int16{1};
    vars.int16.i = vars_int16{2};
    vars.int16.i_next = vars_int16{3};
    vars.int16.id = vars_int16{4};
    vars.int16.id_next = vars_int16{5};

    % Build a copy from the variable data
    obs2_restored = set_obs_vars(MKF_SF1.copy(), vars);  % makes a new copy
    assert(isequal(obs2_restored.xkp1_est, obs2.xkp1_est))
    %assert(isequal(obs2_restored.Pkp1, obs2.Pkp1))
    assert(isequal(obs2_restored.p_seq_g_Yk, obs2.p_seq_g_Yk))
    assert(isequal(obs2_restored.filters.Xkp1_est, obs2.filters.Xkp1_est))
    assert(isequal(obs2_restored.filters.Pkp1, obs2.filters.Pkp1))
    assert(isequal(obs2_restored.rk, obs2.rk))
    assert(isequal(obs2_restored.i, obs2.i))
    assert(isequal(obs2_restored.i_next, obs2.i_next))

    % Observer updates
    obs2.update(yk, uk);

    % Convert dynamic variables to vdata struct
    vars = get_obs_vars(obs2);
    vars_double = {vars.xkp1_est, vars.p_seq_g_Yk, ...
        vars.xkp1_est_f, vars.Pkp1_f};
    vdata = make_data_vectors(vars_double);
    vdata_int16 = make_data_vectors(struct2cell(vars.int16), 'int16');

    % Save observer estimates
    Xk_est{2}(i+1,:) = obs2.xk_est';
    Yk_est{2}(i+1,:) = obs2.yk_est';

end

% Check observer estimates are identical
assert(isequaln(Xk_est{1}, Xk_est{2}))
assert(isequaln(Yk_est{1}, Yk_est{2}))