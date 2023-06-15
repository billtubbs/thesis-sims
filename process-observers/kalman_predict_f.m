function [xkp1_est,Pkp1] = kalman_predict_f(A,B,Q,xk_est,Pk,uk)
% [xkp1_est,Pkp1] = kalman_predict_f(A,B,Q,xk_est,Pk,uk)
% computes the predicted states, outputs, and state error 
% covariance of a Kalman filter (filtering form) in the next
% time instant given the model parameters, the current error
% covariance, and known input:
%
%   x_est(k+1|k) = A(k) * x_est(k|k) + B(k) * u(k)
%   P(k+1|k) = A(k) * P(k|k) * A'(k) + Q(k)
%
% Arguments:
%   A, B : matrices
%       Discrete-time system model matrices (assumes D = 0).
%   Q : matrix, size (n, n)
%       State error covariance.
%   xk_est : (n, 1) vector
%     Posterior estimate of the states, x_est(k|k).
%   Pk : (n, n) matrix
%     State error covariance, P(k|k).
%   uk : (nu, 1) vector
%     Known system inputs at time k.
%
% Returns:
%   xkp1_est : (n, 1) vector
%     Predicted states at next time instant, x_est(k+1|k).
%   Pk : (n, n) matrix
%     Predicted state error covariance at next time instant,
%     P(k+1|k).
%

    % Make predictions of states and outputs in next time step
    % x_est(k+1|k) and y_est(k+1|k) i.e. based on the data up
    % to time k. These will be used as priors in the next time 
    % step.
    xkp1_est = A * xk_est + B * uk;

    % Error covariance of predicted state estimates
    Pkp1 = A * Pk * A' + Q;

end