function params = get_obs_params(obs)
% T = get_obs_params(obs) returns a struct containing
% selected parameters of the observer. Which params
% are selected depends on the observer label.

    if strcmp(obs.label, 'none')  % no observer

        params = struct();

    elseif startsWith(obs.label, 'KF')  % Standard Kalman filters

        % Params to return
        params.P0 = obs.P0;
        params.Q = obs.Q;
        params.R = obs.R;

    elseif startsWith(obs.label, 'SKF')  % Scheduled Kalman filters

        % Params to return
        params.P0 = obs.P0;
        params.R = obs.R;
        params.Q0 = obs.Q0;
        params.sigma_wp = obs.sigma_wp;

    elseif startsWith(obs.label, 'MMKF')  % adaptive multi-model Kalman filter

        % Params to return
        params.P0 = obs.P0;
        params.Q0 = obs.Q0;
        params.R = obs.R;
        params.epsilon = obs.epsilon;
        params.sigma_wp = obs.sigma_wp;
        params.n_filt = obs.n_filt;
        params.f = obs.f;
        params.n_min = obs.n_min;

    else
        error('Value error: observer type not recognized')
    end

end