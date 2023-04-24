function sim_out = run_obs_simulation(Ts, data, observers)
% sim_out = run_obs_simulation(Ts, data, observers)
% Runs a simulation of a set of pre-defined observers fed with
% input data U and output data Y.
%
% Arguments:
%   Ts : sample time.
%   data : table of data including measured input signals named
%       'U', measured output signals named 'Y'. May also contain
%       data for other variables to be used by the observers
%       or to be inlcuded in the simulation outputs.
%   observers : cell array of observer structs.
%
% Returns:
%   sim_out : struct with the following fields: 
%       data : table containing simulation results.
%       if one or more multi-model filters are included:
%       MKF_obs : array of indices indicating which 
%           observers are multi-model observers.
%       MKF_i : filter sequence index.
%       MKF_X_est : array of estimates of each Kalman filter.
%       MKF_p_seq_g_Yk : conditional probabilities.
%

    % Extract data needed for simulation
    nT = size(data, 1) - 1;
    U = data{:, 'U'};
    Y_m = data{:, 'Y_m'};

    % Labels
    data_labels = data.Properties.VariableNames;
    other_data_labels = setdiff(data_labels, {'U', 'Y_m'});

    % Simulation setup
    k = (0:nT)';
    t = Ts*k;

    % Prepare arrays for storing simulation variables
    n_obs = numel(observers);
    X_est = cell(1, n_obs);
    Y_est = cell(1, n_obs);
    for i = 1:n_obs
        n = observers{i}.n;
        ny = observers{i}.ny;
        X_est{i} = nan(nT+1, n);
        Y_est{i} = nan(nT+1, ny);
    end

    % Find and count number of MMKF filters
    n_obs_mkf = 0;
    obs_mkf = nan(1, n_obs_mkf);
    for i = 1:n_obs
        if ismember(observers{i}.type,  ...
                ["MKF" "MKF_SF_RODD95" "MKF_SF_RODD" "MKF_SP_RODD"] ...
            )
            n_obs_mkf = n_obs_mkf + 1;
            obs_mkf(n_obs_mkf) = i;
        end
    end

    % Cell arrays to store data on each MKF observer
    MKF_i = cell(1, n_obs_mkf);
    MKF_X_est = cell(1, n_obs_mkf);
    MKF_Y_est = cell(1, n_obs_mkf);
    MKF_p_seq_g_Yk = cell(1, n_obs_mkf);

    for f = 1:n_obs_mkf
        obs = observers{obs_mkf(f)};
        MKF_i{f} = int16(zeros(size(k)));
        switch obs.type
            case {"MKF", "MKF_SP", "MKF_SP_RODD"}
                nh = obs.nh;
            case {"MKF_SF_RODD95", "MKF_SF_RODD"}
                nh = obs.nm;  % save merged estimates
        end
        MKF_X_est{f} = nan(size(k, 1), nh*n);
        MKF_p_seq_g_Yk{f} = nan(size(k, 1), nh);
    end

    % Start simulation at k = 0 (i = 1)
    for ki = 0:nT
        
        i = ki + 1;  % array index starts at 1
        
        % Get current process inputs
        uk_m = U(i, :)';

        % Get current process outputs
        yk_m = Y_m(i, :)';

        % Update observers using measurements uk_m and yk_m
        for j = 1:n_obs
            obs = observers{j};

            switch obs.type
                case "KFF"  % Kalman filter - filtering form

                    % Update observer estimates
                    obs.update(yk_m, uk_m);

                case "SKF"  % Scheduled Kalman filter

                    % Provide actual shock occurence
                    rk = data{i, 'gamma'} + 1;  % 1-based indexing
    
                    % Update observer estimates
                    obs.update(yk_m, uk_m, rk);
                
                case {"MKF", "MKF_SP", "MKF_SF_RODD95", "MKF_SF_RODD", ...
                        "MKF_SP_RODD"}
                    % Multi-model observers with pruning algorithm

                    % Update observer estimates
                    obs.update(yk_m, uk_m);
    
                    % Save additional simulation data for plotting
                    f_mkf = find(obs_mkf == j);
                    if isprop(obs, "i")
                        MKF_i{f_mkf}(i) = obs.i(1);
                    else
                        MKF_i{f_mkf}(i) = nan;
                    end
                    switch obs.type
                        case {"MKF", "MKF_SP", "MKF_SP_RODD"}
                            MKF_X_est{f_mkf}(i, :) = ...
                                reshape(obs.filters.Xk_est, 1, []);
                            MKF_Y_est{f_mkf}(i, :) = ...
                                reshape(obs.filters.Yk_est, 1, []);
                            MKF_p_seq_g_Yk{f_mkf}(i, :) = ...
                                obs.p_seq_g_Yk';
                        case {"MKF_SF_RODD95", "MKF_SF_RODD"}
                            % For these observers save merged 
                            % estimates
                            MKF_X_est{f_mkf}(i, :) = ...
                                reshape(obs.merged.Xk_est, 1, []);
                            MKF_Y_est{f_mkf}(i, :) = ...
                                reshape(obs.merged.Yk_est, 1, []);
                            MKF_p_seq_g_Yk{f_mkf}(i, :) = ...
                                obs.merged.p_seq_g_Yk';
                    end

                otherwise
                    error('Value error: observer type not recognized')

            end

        end

        % Record observer estimates
        for j = 1:n_obs
            X_est{j}(i, :) = observers{j}.xk_est';
            Y_est{j}(i, :) = observers{j}.yk_est';
        end

    end

    % Calculate output estimation errors
    E_obs = repmat(Y_m, 1, n_obs) - cell2mat(Y_est);

    % Convert estimates to tables
    X_est = table(X_est{:}, 'VariableNames', ...
        vector_element_labels('X_est', '', n_obs));
    Y_est = table(Y_est{:}, 'VariableNames', ...
        vector_element_labels('Y_est', '', n_obs));

    % Main simulation results table
    sim_out.data = [table(k,t,U,Y_m) X_est Y_est table(E_obs) ...
        data(:, other_data_labels)];

    % Data on MKF filters
    sim_out.MKF_obs = obs_mkf;
    sim_out.MKF_i = MKF_i;
    sim_out.MKF_X_est = MKF_X_est;
    sim_out.MKF_Y_est = MKF_Y_est;
    sim_out.MKF_p_seq_g_Yk = MKF_p_seq_g_Yk;
    sim_out.observers = observers;

end