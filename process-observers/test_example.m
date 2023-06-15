% Test Minimal example used in ReadMe

clear all

% Known inputs
U = [     0     0     1     1     1     1     1     1     1     1 ...
          1     1     1     1     1     1     1     1     1     1 ...
          1]';

% Output measurements
Ym =  [    0.2688    0.9169   -1.1294    0.7311    0.6694 ...
           0.0032    0.5431    1.0032    2.6715    2.3024 ...
           0.2674    2.4771    1.3345    0.9487    1.3435 ...
           0.8878    0.9311    1.7401    1.7012    1.7063 ...
           1.3341]';

% Sampling period
Ts = 0.5;

% Discrete-time transfer function
Gpd = tf(0.3, [1 -0.7], Ts);

% State-space representation of above process model
model.A = 0.7;
model.B = 1;
model.C = 0.3;
model.Ts = Ts;
[n, nu, ny] = check_model(model);

% Kalman filter parameters
P0 = 1;  % estimated variance of the initial state estimate
model.Q = 0.01;  % estimated process noise variance
model.R = 0.5^2;  % estimated measurement noise variance

% Kalman filter 1 - prediction form
KF1 = KalmanFilterP(model,P0,'KF1');

% Kalman filter 2 - filtering form
KF2 = KalmanFilterF(model,P0,'KF2');


%% Simulate the observer and record the output estimates:

% Number of sample periods
nT = size(Ym, 1) - 1;
% Arrays to store observer estimates
Xk_est = nan(nT+1, n);
Yk_est = nan(nT+1, 1);
Xkp1_est = nan(nT+1, n);
Ykp1_est = nan(nT+1, 1);

% Save initial estimate (at t=0)
Xkp1_est(1,:) = KF1.xkp1_est';
Ykp1_est(1,:) = KF1.ykp1_est;
for i = 1:nT

    % Update observers with measurements
    KF1.update(Ym(i), U(i));
    KF2.update(Ym(i), U(i));

    % Prediction of states and outputs at next sample time
    Xkp1_est(i+1,:) = KF1.xkp1_est';
    Ykp1_est(i+1,:) = KF1.ykp1_est;

    % Check preedictions are identical
    assert(abs(KF1.xkp1_est - KF2.xkp1_est) < 1e-14)
    assert(abs(KF1.Pkp1 - KF2.Pkp1) < 1e-14)

    % Estimate of states and output at current time
    Xk_est(i,:) = KF2.xk_est;
    Yk_est(i,:) = KF2.yk_est;

end

% Check results

Ykp1_est_test = [
         0    0.0498    0.1063    0.3245    0.5359    0.6769    0.7679 ...
    0.8360    0.8862    0.9298    0.9578    0.9671    0.9844    0.9908 ...
    0.9933    0.9970    0.9974    0.9979    1.0021    1.0049    1.0068 ...
]';
assert(isequal(round(Ykp1_est, 4), Ykp1_est_test))

Yk_est_test = [
    0.0712    0.1518    0.0350    0.3370    0.5384    0.6685    0.7658 ...
    0.8374    0.8997    0.9398    0.9529    0.9777    0.9868    0.9905 ...
    0.9958    0.9963    0.9969    1.0030    1.0070    1.0098       NaN
]';
assert(isequaln(round(Yk_est, 4), Yk_est_test))

% Plot observer output estimates to measurement data

% figure(1)
% t = Ts*(0:nT)';
% plot(t,Ym,'o',t,Ykp1_est,'.-',t,Yk_est,'.-')
% grid on
% xlabel('Time')
% ylabel('Process output')
% legend('y(k)','y(k+1) prediction','y(k) estimate')
% title("Observer estimates compared to process measurements")