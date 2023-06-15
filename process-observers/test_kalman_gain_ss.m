% Test static Kalman filter gain calculation
%
%  - kalman_gain_ss.m
%

clear all

%% Example from course notes, exercise 12.3

% System model
A = 0.8; 
B = 2; 
C = 0.5; 
Q = 0.3; 
R = 1.1;

[Kf, P] = kalman_gain_ss(A, C, Q, R);

% Check P solves the Riccati equation
assert(abs(A*(P - P*C'*pinv(C*P*C' + R)*C*P)*A' + Q - P) < 1e-15)

% Gain for prediction form of KF
K = A * Kf;
assert(round(P, 4) == 0.6741)
assert(round(K, 4) == 0.2126)

K_calc = A*P*C'*(C*P*C' + R)^-1;
assert(abs(K_calc - K) < 1e-16)

% Gain for filtering form of KF
Kf_calc = P*C'*(C*P*C' + R)^-1;
assert(Kf_calc == Kf)


%% Example system adapted from MATLAB documentation

% System model
A = [-0.9, 0.7; 
     -0.3, 0.1];
C = [1, 1];
Q = [1, 0;
     0, 3];
R = 0.1;

[Kf, P] = kalman_gain_ss(A, C, Q, R);
K = A * Kf;

assert(isequal(round(P, 4), [ ... 
    4.7687    0.9438
    0.9438    3.2369 ...
]))
assert(isequal(round(K, 4), [ ... 
   -0.2216; ...
   -0.1297 ...
]))

K_calc = A * P* C' * (C * P * C' + R)^-1;
assert(all(abs(K_calc - K) < 1e-15))

% Filtering form
Kf_calc = P * C' * pinv(C * P * C' + R);
assert(isequal(round(Kf_calc, 4), [ ... 
   0.5716; ...
   0.4184 ...
]))
assert(all(abs(Kf_calc - Kf) < 1e-15))

% Using the dlqe function
G = eye(2);  % apply process noises to all states
[Kf_calc2,~,~] = dlqe(A,G,C,Q,R);
assert(all(abs(Kf_calc - Kf_calc2) < 1e-15))

%% Example system from Mathworks answer

% System model
A = [0.7 1;
     0   1];
B = [1; 0];
C = [0.3 0];
D = 0;

% Error covariance matrices
Q = diag([0.0001 0.01]);
R = 0.01;

[Kf, P] = kalman_gain_ss(A, C, Q, R);
K = A * Kf;

assert(isequal(round(P, 4), [ ... 
    0.0857    0.0444
    0.0444    0.0368 ...
]))
assert(isequal(round(Kf, 4), [ ... 
   1.4515; ...
   0.7514 ...
]))

% Check results satisfy eqn.s
assert(all(A*(P - P*C'*(C*P*C' + R)^-1*C*P)*A' + Q - P < 1e-14, [1 2]))
assert(all(K - A*P*C'*(C*P*C' + R)^-1 < 1e-14, [1 2]))

% Using MATLAB Kalman function

% Construct model
n = 2; ny = 1; Ts = 1;
N = zeros(n, ny);
G = eye(n);  % apply process noises to all states
H = zeros(ny, n);  % no direct transmission of noises
Gmodel = ss(A, [B G], C, [D H], Ts);

% Use MATLAB's Kalman filter function
[~, K2, P2] = kalman(Gmodel, Q, R, N, 'delayed');
assert(all(K - K2 < 1e-14, [1 2]))
assert(all(P - P2 < 1e-14, [1 2]))
