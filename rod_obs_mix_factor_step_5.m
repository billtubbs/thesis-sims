%% System model definition
%
% Discrete transfer function polynomial model identified
% from SAG mill circuit simulation model.
% 
% Usage: run this script from a main script to define
%  - Gd : Process transfer function Gd(z)
%  - HNd : ARIMA noise disturbance transfer function HNd(z)
%  - HDd : Input disturbance RODD transfer function
%  - Gpd : combined system transfer function (2 inputs, 1 output)
%  - Gpss : state space model of combined system
%

%% Discrete transfer function polynomial models

% Sample time (hrs)
Ts = 3/60;

% Process model (from system identification)
% P1
%Gc = tf(-34.3, [0.429 1], 'TimeUnit', 'minutes');
% P2U
%Gc = tf(-32.30, [0.06961 0.3407 1], 'TimeUnit', 'minutes');
% P1D_c4
s = tf('s');
Gc = -33.4 * exp(-0.15 * s) / (1 + 0.273*s);
Gc.TimeUnit = 'minutes';
Gd = c2d(Gc, Ts, 'ZOH');
Gdss = ss(absorbDelay(Gd));

% RODD step disturbance process
ThetaD = 1;
PhiD = 1;
d = 1;
HDd = rodd_tf(ThetaD, PhiD, d, Ts);
HDd.TimeUnit = 'minutes';

% Combined system transfer functions
Gpd = series(Gd, HDd);

% State space representation
Gpss2 = minreal(absorbDelay(ss(Gpd)));
% A = Gpss.A;
% B = Gpss.B;
% C = Gpss.C;
% D = Gpss.D;

% Discrete time state space model - P1D_c4
A = [  0.8326       1       0       0       0;
            0       0       1       0       0;
            0       0       0       1       0;
            0       0       0       0       1;
            0       0       0       0       1];
B = [   0
        0
        0
        0
        1];
C = [    -5.59       0       0       0       0];
D = 0;
Gpss = ss(A, B, C, D, Ts, 'TimeUnit', 'minutes');

% % Discrete time state space model - P2U
% A = [   2.75116  -1.26704  0.782923
%              2        0        0
%              0      0.5        0];
% B = [1;
%      0;
%      0];
% C = [0  -0.266944   -0.492035];
% D = 0;
% Gpss = ss(A, B, C, D, Ts, 'TimeUnit', 'minutes');

% Dimensions
n = size(A, 1);
nu = size(B, 2);
ny = size(C, 1);

% Check both versions are almost identical
wp = idinput(15, 'PRBS');
t_test_sim = Ts*(0:14)';
[y1, ~] = lsim(Gpd, wp, t_test_sim);
[y2, ~] = lsim(Gpss, wp, t_test_sim);
[y3, ~] = lsim(Gpss2, wp, t_test_sim);
assert(max(abs(y1 - y2)) < 0.05)
assert(max(abs(y1 - y3)) < 0.05)

% Designate which input and output variables are
% measured
u_meas = false;
y_meas = true;

% Default initial condition
x0 = zeros(n, 1);


%% Parameters for random inputs

% RODD random variable parameters
epsilon = 0.01;
sigma_wp = [0.002717 0.2717];  % [0.01 1] .* step magnitude

% Process noise standard deviation
sigma_W = [0; 0];

% Measurement noise standard deviation
sigma_M = 5;

% To check observer with no noise disturbances
%sigma_W = [0; 0];
%sigma_M = 0.001;

% Initial state of disturbance process
p0 = 0;


%% Define observers

% Observer model (without disturbance noise input)
Bu = B(:, u_meas);
Du = D(:, u_meas);
nu = sum(u_meas);

% Disturbance input (used by SKF observer)
Bw = B(:, ~u_meas);
nw = sum(~u_meas);

% Check observability of system
Qobs = obsv(A,C);
unobs = length(A) - rank(Qobs);
assert(unobs == 0);

% Common parameters for all observers
P0 = diag([0.1 0.1 0.1 0.1 0.01]);  % initial error covariance
R = sigma_M^2;

% process noises for state x1
q00 = 0.01^2;

% Kalman filter 1 - tuned to sigma_wp(1)
% Covariance matrices
Q = diag([q00 q00 q00 q00 sigma_wp(1)^2]);
KF1 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF1');

% Kalman filter 2 - tuned to sigma_wp(2)
% Covariance matrices
Q = diag([q00 q00 q00 q00 sigma_wp(2)^2]);
KF2 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF2');

% Kalman filter 3 - manually tuned
% Covariance matrices
Q = diag([q00 q00 q00 q00 0.027^2]);
KF3 = kalman_filter(A,Bu,C,Du,Ts,P0,Q,R,'KF3');

% Scheduled Kalman filter
% Use this to test performence of multi-model filters
Q0 = diag([q00 q00 q00 q00 0]);
SKF = kalman_filter(A,Bu,C,Du,Ts,P0,Q0,R,'SKF');
SKF.Q0 = Q0;
SKF.Bw = Bw;
SKF.sigma_wp = sigma_wp;

% Multiple model AFMM filter 1
label = 'AFMM1';
Q0 = diag([q00 q00 q00 q00 0]);
f = 100;  % sequence history length
n_filt = 10;  % number of filters
n_min = 5;  % minimum life of cloned filters
AFMM1 = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);

% Multiple model AFMM filter 2
label = 'AFMM';
Q0 = diag([q00 q00 q00 q00 0]);
f = 100;  % sequence history length
n_filt = 20;  % number of filters
n_min = 10;  % minimum life of cloned filters
AFMM = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
    Q0,R,n_filt,f,n_min,label);
