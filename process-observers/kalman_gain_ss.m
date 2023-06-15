function [Kf, P] = kalman_gain_ss(A, C, Q, R)
% [Kf, P] = kalman_gain_ss(A, C, Q, R)
% Compute the steady-state correction gain and error 
% covariance matrix for the filtering form of the 
% Kalman Filter:
%
%   Kf = P * C' * (C * P * C' + R)
%
% Solves the algebraic Riccati equation:
%
%   P = A*(P - P*C'*(C*P*C' + R)^-1*C*P)*A' + Q
%
% For the prediction form of the Kalman filter, the
% correction gain is
%
%   K = A * Kf;
%

    % Using the Matlab idare function:
    %[P, K, ~, ~] = idare(A', C', Q, R, [], []);
    %K = K';
    %Kf = A \ K;

    % Using dlqe function from the control systems toolbox:
    % Note: this function is marked as 'obsolete'
    n = size(A, 1);
    [Kf, P, ~] = dlqe(A, eye(n), C, Q, R);

end
