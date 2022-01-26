function obs = update_SKF(obs, uk, yk, alpha_k)
% obs = update_SKF(obs, uk, yk, alpha)
% updates the scheduled multi-model Kalman filter and 
% calculate the estimates of the states and output at the 
% next sample time. Note: this is a hypothetical observer
% intended only for benchmarking purposes.
%
% Arguments:
%   obs : struct containing the multi-model Kalman filter
%       variables (see function mkf_filter).
%   uk : vector (nu, 1) of system inputs at the current 
%       sample time.
%   yk : vector (ny, 1) of system output measurements
%       at the current sample time.
%   alpha_k : row vector indicating the actual random 
%       shock occurence at the current sample tim.
%

    % Set process noise covariance matrix Q based on actual 
    % shock occurence
    n_dist = size(alpha_k, 2);
    x_var = diag(obs.Q0);
    w_var = (obs.sigma_wp(sub2ind(size(obs.sigma_wp), 1:n_dist, ...
        alpha_k+1)).^2)';
    x_var = x_var + obs.Bw * w_var;
    obs.Q = diag(x_var);

    % Update observer estimates
    obs = update_KF(obs, uk, yk);

end