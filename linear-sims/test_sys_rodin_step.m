%% Test sys_rodin_step.m

clear all
rng(0)

sys_rodin_step

% Check model parameters haven't changed
assert(isequal(Ts, 0.5));
assert(isequal(A, [0.7 1; 0 1]));
assert(isequal(B, [1 0; 0 1]));
assert(isequal(C, [0.3 0]));
assert(isequal(D, zeros(1, 2)));

% Check models are identical
t_test = 0.5*(0:10)';
U_test = reshape(idinput(11*2), 11, 2);
[y_test, t_test] = lsim(Gpss, U_test, t_test);
[y_test2, t_test] = lsim(Gpd, U_test, t_test);
assert(all(abs(y_test - y_test2) < 1e-6))

assert(isequal(u_known, [true; false]));
assert(isequal(y_meas, true));

% Check disturbance parameters
assert(epsilon == 0.01);
assert(isequal(sigma_wp, {[0.01 1]}));

% Process noise standard deviation
assert(isequal(sigma_W, [0; 0]));

% Measurement noise standard deviation
assert(sigma_M == 0.1);

assert(isequal(x0, [0; 0]));
assert(p0 == 0);