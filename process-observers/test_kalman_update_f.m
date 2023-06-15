% Test alternative Kalman filter update equations


%% 2x2 example system with random values

% Discrete time state space model
A = [ 0.8890       0     1 -0.5;
           0  0.8890  -0.5    1;
           0       0     1    0;
           0       0     0    1];
B = [    1 -0.5  0  0;
      -0.5    1  0  0;
         0    0  1  0;
         0    0  0  1];
C = [ 0.1110       0  0  0;
           0  0.1110  0  0];
D = zeros(2, 4);
n = 4;

rng(0);
xk_pred = randn(n, 1);
Pk_pred = randn(n).^2;
sigma_M = [0.5 0.3]';
R = diag(sigma_M.^2);
yk = C * xk_pred + sigma_M .* randn(2, 1);

% Calculation method I was using initially

% Error covariance of output prediction error
Sk = C * Pk_pred * C' + R;

% Correction gain
Kf = Pk_pred * C' / Sk;

% Update state predictions (prior estimates) using 
% measurements from current time to produce 'a posteriori' 
% state estimates
xk_est = xk_pred + Kf * (yk - C * xk_pred);

% Updated output estimate
yk_est = C * xk_est;

% Updated error covariance of state estimates
Pk1 = (eye(size(Pk_pred)) - Kf * C) * Pk_pred;

% Alternative calculation method
% The following is the 'Joseph Stabilized' version of
% the update equation which guarantees that Pk is positive
% semi-definite in the presence of roundoff error.
% (See p73 of Lewis etal. Book Optimal and Robust Estimation).
z = eye(size(Pk_pred)) - Kf * C;
Pk = z * Pk_pred * z' + Kf * R * Kf';

% Make sure difference is small
assert(max(abs(Pk - Pk1), [], [1 2]) < 1e-14)


%% 4x4 system with values from actual simulation

% When I ran a simulation using the initial method
% of calculating Pk, I found that the divergence between
% this method and the stabilized method increased over
% the course of the simulation. In this example, I took
% values of xk_pred, Pk_pred from the simulation when
% the maximum difference between the two calculations was
% greater than 1e-4.

load("data/sim_vars_for_testing.mat")

% Error covariance of output prediction error
Sk = C * Pk_pred * C' + R;

% Correction gain
Kf = Pk_pred * C' / Sk;

% Update state predictions (prior estimates) using 
% measurements from current time to produce 'a posteriori' 
% state estimates
xk_est = xk_pred + Kf * (yk - C * xk_pred);

% Updated output estimate
yk_est = C * xk_est;

% Updated error covariance of state estimates
Pk1 = Pk_pred - Kf * C * Pk_pred;
Pk2 = Pk_pred - (Pk_pred * C' / (C * Pk_pred * C' + R)) * C * Pk_pred;
assert(isequal(Pk1, Pk2))

% Alternative calculation method
% The following is the 'Joseph Stabilized' version of
% the update equation which guarantees that Pk is positive
% semi-definite in the presence of roundoff error.
% (See p73 of Lewis etal. Book Optimal and Robust Estimation).
z = eye(size(Pk_pred)) - Kf * C;
Pk = z * Pk_pred * z' + Kf * R * Kf';



