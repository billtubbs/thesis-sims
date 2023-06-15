% Test construct_Q_model_SP.m function
% (used by MKFObserverSP)

clear all


%% Test on SISO system with no RODDs
% with one known input and one standard (Gaussian) input disturbance

Q0 = [ 0.0100         0
            0         0];
B = [1 0; 0 1];
u_known = [true false]';
alpha = 0.01;
sigma_wp = {0.01};

[Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp);
Q_test = {
    [ 0.0100       0
           0       0.0001]
};

assert(isequal(Q, Q_test))
assert(isequal(p_rk, 1))



%% Test on SISO system with one RODD
% with one known input and one randomly-occurring step input disturbance

Q0 = [ 0.0100         0
            0         0];
B = [1 0; 0 1];
u_known = [true false]';
alpha = 0.01;
sigma_wp = {[0.01 1]};

[Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp);
Q_test = {
    [ 0.0100       0
           0       0.0001] ... 
    [ 0.0100       0
           0       1     ]
};

assert(isequal(Q, Q_test))
assert(isequal(p_rk, [0.99; 0.01]))


%% Test on 2x2 system with one RODD

Q0 = [
    0.0100         0         0         0
         0    0.0100         0         0
         0         0         0         0
         0         0         0         0
];

B = [
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
];
u_known = [true true false false]';
epsilon = [0.01];
sigma_wp = {0.01, [0.01 1]};
d = 1;
alpha = (1 - (1 - epsilon).^d);

[Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp);
Q_test = {
    [    0.0100         0         0         0
              0    0.0100         0         0
              0         0    0.0001         0
              0         0         0    0.0001] ...
    [    0.0100         0         0         0
              0    0.0100         0         0
              0         0    0.0001         0
              0         0         0    1.0000]
};
assert(isequal(Q, Q_test));
assert(isequal(round(p_rk, 6), [0.99 0.01]'));

% Test with non-zero covariances of process states
Q0 = [
    0.0105     0.002         0         0
     0.003    0.0090         0         0
         0         0         0         0
         0         0         0         0
];
[Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp);
Q_test = {
    [    0.0105     0.002         0         0
          0.003    0.0090         0         0
              0         0    0.0001         0
              0         0         0    0.0001] ...
    [    0.0105     0.002         0         0
          0.003    0.0090         0         0
              0         0    0.0001         0
              0         0         0    1.0000]
};
assert(isequal(Q, Q_test));
assert(isequal(round(p_rk, 6), [0.99 0.01]'));


%% Test on 2x2 system with two RODDs

Q0 = [
    0.0100         0         0         0
         0    0.0100         0         0
         0         0         0         0
         0         0         0         0
];

B = [
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
];
u_known = [true true false false]';
epsilon = [0.01; 0.01];
sigma_wp = {[0.01 1], [0.01 1]};
d = 1;
alpha = (1 - (1 - epsilon).^d);

[Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp);
Q_test = {
    [    0.0100         0         0         0
              0    0.0100         0         0
              0         0    0.0001         0
              0         0         0    0.0001] ...
    [    0.0100         0         0         0
              0    0.0100         0         0
              0         0    1.0000         0
              0         0         0    0.0001] ...
    [    0.0100         0         0         0
              0    0.0100         0         0
              0         0    0.0001         0
              0         0         0    1.0000]
};
assert(isequal(Q, Q_test));
assert(isequal(round(p_rk, 6), [0.9801 0.0099 0.0099]'));


%% Test on 3x2 system with no RODDs

Q0 = [
    0.0100         0         0         0         0
         0    0.0100         0         0         0
         0         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0
];

B = [
     1     0     0     0     0
     0     1     0     0     0
     0     0     1     0     0
     0     0     0     1     0
     0     0     0     0     1
];
u_known = [true true false false false]';
epsilon = [];
sigma_wp = {0.01, 0.02, 0.03};
d = 1;
alpha = (1 - (1 - epsilon).^d);

[Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp);
Q_test = {
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0    0.0001         0         0
              0         0         0    0.0004         0
              0         0         0         0    0.0009] ...
};
assert(isequal(Q, Q_test));
assert(p_rk == 1)


%% Test on 3x2 system with 3 RODDs

Q0 = [
    0.0100         0         0         0         0
         0    0.0100         0         0         0
         0         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0
];

B = [
     1     0     0     0     0
     0     1     0     0     0
     0     0     1     0     0
     0     0     0     1     0
     0     0     0     0     1
];
u_known = [true true false false false]';
epsilon = [0.005; 0.0025; 0.0025];
sigma_wp = {
    [0.0100    1.0000],
    [0.0050    0.5000],
    [0.0050    0.5000]
};
d = 1;
alpha = (1 - (1 - epsilon).^d);

[Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp);
Q_test = {
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0    0.0001         0         0
              0         0         0    2.5e-5         0
              0         0         0         0    2.5e-5] ...
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0         1         0         0
              0         0         0    2.5e-5         0
              0         0         0         0    2.5e-5] ...
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0    0.0001         0         0
              0         0         0      0.25         0
              0         0         0         0    2.5e-5] ...
    [    0.0100         0         0         0         0
              0    0.0100         0         0         0
              0         0    0.0001         0         0
              0         0         0    2.5e-5         0
              0         0         0         0      0.25]
};
assert(isequal(Q, Q_test));

p_rk_test1 = [
    prod(1-epsilon);
    epsilon(1)*prod(1-epsilon(2:3));
    epsilon(2)*prod(1-epsilon([1, 3]));
    epsilon(3)*prod(1-epsilon(1:2));
];
p_rk_test2 = [0.990031 0.004975 0.002481 0.002481]';
assert(max(abs(p_rk - p_rk_test1)) < 1e-12)
assert(isequal(round(p_rk, 6), p_rk_test2))
