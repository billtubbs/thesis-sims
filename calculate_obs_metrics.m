function [metrics, params, errors, labels] = calculate_obs_metrics(Y, ...
        Y_est, obs_labels, Pd, Ts, tau_ss)
% Used by rod_obs_sim.m to compute observer metrics from
% simulation results.
%

    % Labels to be used
    labels = {'RMSE', 'RMSE_tr', 'RMSE_ss', 'Var_ss', 'RMSD_ss'};

    % Store results in struct
    params = struct();
    errors = struct();

    % Determine number of results
    n_obs = length(obs_labels);
    assert(rem(size(Y_est, 2), n_obs) == 0)
    % TODO: currently does not support ny > 1
    %ny = size(Y_est, 2) / ny;

    % Errors in observer state estimates
    errors.Y = repmat(Y, 1, n_obs) - Y_est;

    % Mean-squared errors
    Y_RMSE = sqrt(mean(errors.Y.^2, 1));

    % Calculate metrics for steady-state periods and after steps
    params.tau_ss = tau_ss;  % settling time
    params.n_settle = ceil(params.tau_ss/Ts);
    errors.ss_periods = steady_state_periods(Pd, params.n_settle);
    Y_RMSE_tr = sqrt( ...
        mean(errors.Y(~errors.ss_periods,:).^2, 1) ...
    );
    Y_RMSE_ss = sqrt( ...
        mean(errors.Y(errors.ss_periods,:).^2, 1) ...
    );

    % Record number of samples in MSE calculations
    params.nT_Y_RMSE = size(errors.Y, 1);
    params.nT_Y_RMSE_tr = sum(~errors.ss_periods);
    params.nT_Y_RMSE_ss = sum(errors.ss_periods);
    assert(params.nT_Y_RMSE_tr + params.nT_Y_RMSE_ss == params.nT_Y_RMSE)

    % Calculate variance of estimates during steady-state periods
    trans_idxs = transition_periods(Pd);
    n_resp = numel(trans_idxs);
    sq_devs = nan(size(Y_est));
    for i = 1:n_resp
        idx = trans_idxs{i};
        idx1 = idx(1) + params.n_settle - 1;
        avgs = mean(Y_est(idx1:idx(2), :));
        sq_devs(idx1:idx(2), :) = (Y_est(idx1:idx(2), :) - avgs).^2;
    end
    sq_devs = sq_devs(~isnan(sq_devs(:, 1)), :);
    Y_var_ss = mean(sq_devs);

    % Calculate mean-squared differences
    diffs_Y = [nan(1, size(Y_est, 2)); diff(Y_est, 1)];
    Y_RMSD_ss = sqrt( ...
        mean(diffs_Y(errors.ss_periods,:).^2, 1, 'omitnan') ...
    );

    % Make summary table
    metrics = array2table( ...
        [Y_RMSE' Y_RMSE_tr' Y_RMSE_ss' Y_var_ss' Y_RMSD_ss'], ...
        'RowNames', string(obs_labels), ...
        'VariableNames', labels ...
    );

end