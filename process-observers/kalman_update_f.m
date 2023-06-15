function [xk_est,Pk,yk_est,Kf,Sk] = kalman_update_f(C,R,xk_pred,Pk_pred,yk)
% [xk_est,Pk,yk_est,Kf,Sk] = kalman_update_f(C,R,xk_pred,Pk_pred,yk)
% computes the updated correction gain, state estimate, state error 
% covariance, output estimate, and output error covariance in the
% current time instant of a Kalman filter (filtering form) given
% a new measurement:
%
%   S(k) = C(k) * P(k|k-1) * C'(k) + R(k)
%   Kf(k) = P(k)*C(k)'*(C(k)*P(k)*C'(k) + R(k))^-1
%   x_est(k|k) = x_est(k|k-1) + Kf * (y(k) - C(k) * x_est(k|k-1))
%   P(k|k) = P(k) - Kf(k)*S(k)*Kf'(k) + Q(k)
%   y_est(k|k) = C(k) * x_est(k|k)
%
% Arguments:
%   C : (ny, n) matrix
%     Measurement matrix.
%   R : (ny, ny) matrix
%     Measurement error covariance.
%   xk_pred : (n, 1) vector
%     Prior estimate of the states, x_est(k|k-1).
%   Pk_pred : (n, n) matrix
%     Prior estimate of the state error covariance, P(k|k-1).
%   yk : (ny, 1) vector
%     System output measurement at time k.
%
% Returns:
%   xk_est : (n, 1) vector
%     Posterior estimate of the states, x_est(k|k).
%   Pk : (n, n) matrix
%     Covariance of state estimation error, P(k|k).
%   yk_est : (ny, 1) vector
%     Posterior estimate of the outputs, y_est(k|k).
%   Kf : (n, ny) matrix
%     Correction gain matrix.
%   Sk : (ny, ny) matrix
%     Covariance of output prediction error, S(k).
%

    % Error covariance of output prediction error
    Sk = C * Pk_pred * C' + R;

    % Correction gain
    Kf = Pk_pred * C' / Sk;

    % Update state predictions (prior estimates) using 
    % measurements from current time to produce 'a posteriori' 
    % state estimates
    xk_est = xk_pred + Kf * (yk - C * xk_pred);

    % Update error covariance of state estimates
    %Pk = Pk_pred - Kf * Sk * Kf';  % inaccurate
    %Pk = Pk_pred - Kf * C * Pk_pred;

    % The following is the 'Joseph Stabilized' version of
    % the update equation which guarantees that Pk is positive
    % semi-definite in the presence of roundoff error.
    % (See p73 of Lewis etal. Book Optimal and Robust Estimation).
    z = eye(size(Pk_pred)) - Kf * C;
    Pk = z * Pk_pred * z' + Kf * R * Kf';

    % TODO: Should this check be here or moved to MKF methods
    % where it is apparently needed?
    if ~isscalar(Pk)
        % Ensure Pk is symmetric
        %if ~issymmetric(Pk)
        %    warning("P(k) was not symmetric");
        %end
        Pk = (Pk + Pk.')/2;
        % Check P(k) is also positive definite
        % This is not necessary if above update equation is used
        %[~,err] = cholcov(Pk,0);
        %assert(err == 0)
    end

    % Updated output estimate
    yk_est = C * xk_est;

end