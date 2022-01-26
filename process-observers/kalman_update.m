function [K, P] = kalman_update(P,A,C,Q,R)
% [K, P] = kalman_update(P,A,C,Q,R) computes the updated
% gain and covariance matrix of a Kalman filter.
%
% K = A*P*C'*(C*P*C' + R)^-1;
% P = A*(P - P*C'*(C*P*C' + R)^-1*C*P)*A' + Q;
%

    % Observer gain
    Mx = P * C' / (C * P * C' + R);
    K = A * Mx;

    % Covariance matrix
    P = A * (P - Mx * C * P) * A' + Q;

end