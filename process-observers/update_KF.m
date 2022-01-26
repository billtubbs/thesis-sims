function obs = update_KF(obs, uk, yk)
% obs = update_KF(obs, uk, yk) updates the gain and 
% covariance matrix of the Kalman filter observer and
% calculates the estimates of the states and output at
% the next sample time.
%
% Arguments:
%   obs : struct containing the Kalman filter variables
%         (see function kalman_filter).
%   uk : vector (nu, 1) of system inputs at current sample time k.
%   yk : vector (ny, 1) of system output measurements at current
%       sample time k.
%

% TODO: Consider switching argument order to (obs, yk, uk) for
% consistency with MATLAB's correct(obj,y) and predict(obj)
% functions for extendedKalmanFilter

    if ~obs.static_gain
        % Update observer gain and covariance matrix
        [obs.K, obs.P] = kalman_update(obs.P, obs.A, obs.C, obs.Q, obs.R);
    end

    % Update state and output estimates for next timestep
    obs.xkp1_est = obs.A * obs.xkp1_est + obs.B * uk + ...
        obs.K * (yk - obs.C * obs.xkp1_est);
    obs.ykp1_est = obs.C * obs.xkp1_est;

end