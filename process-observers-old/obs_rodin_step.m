% Various observers designed for the system defined in:
%  - sys_rodin_step.m

% Observer model without disturbance noise input
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);
nw = sum(~u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% Check observability of system
Qobs = obsv(A,C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Specify covariance for state variable 1
% This is used by all observers
Q1 = 0.01;

% Adjustment factor for Q2 (applies to KF1, KF2, KF3, SKF only)
% logspace(-2, 2, 9)
% 0.0100    0.0316    0.1000    0.3162    1.0000    3.1623   10.0000   31.6228  100.0000
% logspace(-1, 1, 9)
% 0.1000    0.1778    0.3162    0.5623    1.0000    1.7783    3.1623    5.6234   10.0000
% logspace(-1, 1, 17)
% 0.5623    0.7499    1.0000    1.3335    1.7783
adj =  1;

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 adj*sigma_wp(1)^2]);
R = sigma_M^2;
KF1 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF1');

% Kalman filter - tuned to sigma_wp(2)
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 adj*sigma_wp(2)^2]);
R = sigma_M^2;
KF2 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF2');

% Kalman filter - manually tuned
% Covariance matrices
P0 = 1000*eye(n);
Q = diag([Q1 adj*0.1^2]);
R = sigma_M^2;
KF3 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF3');

% Scheduled Kalman filter
P0 = 1000*eye(n);
Q0 = diag([Q1 0]);
R = sigma_M^2;
SKF = kalman_filter(A,Bu,C,Du,Ts,P0,Q0,R,'SKF');
SKF.Q0 = Q0;
SKF.Bw = Bw;
SKF.sigma_wp = sigma_wp;

% Multiple model filter 1
label = 'MMKF1';
P0 = 1000*eye(n);
Q0 = diag([Q1 0]);
R = sigma_M^2;
f = 100;  % sequence history length
n_filt = 5;  % number of filters
n_min = 3;  % minimum life of cloned filters
MMKF1 = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Multiple model AFMM filter 2
label = 'MMKF2';
P0 = 1000*eye(n);
Q0 = diag([Q1 0]);
R = sigma_M^2;
f = 100;  % sequence history length
n_filt = 10;  % number of filters
n_min = 4;  % minimum life of cloned filters
MMKF2 = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);
