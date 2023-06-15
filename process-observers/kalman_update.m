function [K, P] = kalman_update(P,A,C,Q,R)
% [K, P] = kalman_update(P,A,C,Q,R) computes the updated
% gain and covariance matrix of a Kalman filter in
% prediction form:
%
% K = A*P*C'*(C*P*C' + R)^-1;
% P = A*(P - P*C'*(C*P*C' + R)^-1*C*P)*A' + Q;
%

    % Observer gain
    Kf = P * C' / (C * P * C' + R);
    K = A * Kf;

    % Covariance matrix
    P = A * (P - Kf * C * P) * A' + Q;

end