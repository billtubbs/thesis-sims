%% Test sys_rodin_step_2x2.m

clear all
rng(0)

sys_rodin_step_2x2sym2

% Check model parameters haven't changed
assert(isequal(Ts, 1));
assert(isequal(A, [ ...
     0.8890         0    1  0.5
          0    0.8890  0.5    1
          0         0    1    0
          0         0    0    1
]));
assert(isequal(B, [ ...
     1  0.5  0  0;
   0.5    1  0  0;
     0    0  1  0;
     0    0  0  1
]));
assert(isequal(C, [ ...
     0.1110       0  0  0;
          0  0.1110  0  0
]));
assert(isequal(D, zeros(2, 4)));

% Check models are identical
t_test = (0:10)';
U_test = reshape(idinput(11*4), 11, 4);
[y_test, t_test] = lsim(Gpss, U_test, t_test);
[y_test2, t_test] = lsim(Gpd, U_test, t_test);
assert(all(abs(y_test - y_test2) < 1e-3, [1 2]))

assert(isequal(u_known, [true; true; false; false]));
assert(isequal(y_meas, [true; true]));

% Check disturbance parameters
assert(isequal(epsilon, [0.01; 0.01]));
assert(isequal(sigma_wp, {[0.01 1], [0.01 1]}));

% Process noise standard deviation
assert(isequal(sigma_W, zeros(4, 1)));

% Measurement noise standard deviation
assert(isequal(sigma_M, [0.2; 0.2]));

assert(isequal(x0, zeros(4, 1)));
assert(isequal(p0, [0; 0]));