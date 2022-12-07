% Test Minimal example used in ReadMe

clear all

% Measured inputs
u = [     0     0     1     1     1     1     1     1     1     1 ...
          1     1     1     1     1     1     1     1     1     1 ...
          1]';

% Output measurements
y_m = [    0.2688    0.9169   -1.1294    0.7311    0.6694 ...
           0.0032    0.5431    1.0032    2.6715    2.3024 ...
           0.2674    2.4771    1.3345    0.9487    1.3435 ...
           0.8878    0.9311    1.7401    1.7012    1.7063 ...
           1.3341]';

% Sampling period
Ts = 0.5;

% Discrete-time transfer function
Gpd = tf(0.3, [1 -0.7], Ts);

% State-space representation of above process model
A = 0.7;
B = 1;
C = 0.3;
D = 0;

% Kalman filter parameters
P0 = 1000;  % estimated variance of the initial state estimate
Q = 0.01;  % estimated process noise variance
R = 0.5^2;  % estimated measurement noise variance
obs = kalman_filter(A,B,C,D,Ts,P0,Q,R,'KF1');


%% Simulate the observer and record the output estimates:

% Number of sample periods
nT = size(y_m,1) - 1;
% Array to store observer estimates
y_est = nan(nT,1);
% Save initial estimate (at t=0)
y_est(1,:) = obs.ykp1_est;
for i = 1:nT

    % update observer
    obs = update_KF(obs, u(i), y_m(i));

    % get estimate of output at next sample time
    y_est(i+1,:) = obs.ykp1_est;

end

% Check results

y_est_test = [
         0    0.1876    0.2997    0.3680    0.5749    0.7048    0.7832 ...
    0.8460    0.8933    0.9359    0.9625    0.9702    0.9866    0.9924 ...
    0.9944    0.9978    0.9979    0.9982    1.0024    1.0051    1.0070 ...
]';
assert(isequal(round(y_est, 4), y_est_test))


% Plot observer output estimates to measurement data

% figure(1)
% t = Ts*(0:nT)';
% plot(t,y_m,'o',t,y_est,'o-')
% grid on
% xlabel('Time')
% ylabel('Process output')
% legend('y_m(k)','y_est(k)')
% title("Observer estimates compared to process measurements")