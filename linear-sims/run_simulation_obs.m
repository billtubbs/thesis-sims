function sim_out = run_simulation_obs(Ts, input_data, state_fcn, ...
    meas_fcn, params, x0, u_known, observers)
% sim_out = run_simulation_obs(Ts, input_data, state_fcn, meas_fcn, ...
%     params, x0, u_meas, observers)
% Runs a simulation of a dynamic system defined by state the
% transition function and measurement function with input data
% U, measurement noise V, and a set of pre-defined observers.
%
% See main script rod_obs_sim2.m for details.
%
% Arguments:
%   Ts : sample time.
%   input_data : table of data including input signal named 'U',
%       measurement noise signal named 'V', and process noise
%       signal named 'W'. May contain data for other variables,
%       to be used by the observers or to be included in the
%       simulation outputs.
%   state_fcn : discrete time state transition function for the
%       system to be simulated, x(k+1) = f(t, x(k), u(k), dt,
%       params).
%   meas_fcn : output function for the system to be simulated,
%       y(k) = h(t, x(k), u(k), dt, params).
%   params : Additional parameters to be passed to state_fcn
%       and meas_fcn (params may be a struct).
%   x0 : Initial values of process states at t=0.
%   u_known : Boolean vector indicating which inputs are measured
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

    % Extract input data needed for simulation
    nT = size(input_data, 1) - 1;
    U = input_data{:, 'U'};
    V = input_data{:, 'V'};
    W = input_data{:, 'W'};
    
    % Labels
    data_labels = input_data.Properties.VariableNames;
    other_data_labels = setdiff(data_labels, {'U', 'V', 'W'});

    % System dimensions
    n = size(x0, 1);
    u0 = U(1, :)';
    nu = size(u0, 1);
    y0 = meas_fcn(x0, u0, Ts, params);
    ny = size(y0, 1);

    % Check input data
    assert(size(u_known, 1) == nu)
    assert(size(V, 2) == ny)
    assert(size(W, 2) == n)

    % Simulation setup
    k = (0:nT)';
    t = Ts*k;

    % Prepare arrays for storing simulation variables
    X = zeros(length(k), n);
    Y = zeros(length(k), ny);
    Y_m = zeros(size(Y));
    n_obs = numel(observers);
    X_est = nan(size(k,1), n_obs*n);
    P = nan(size(k,1), n_obs*n);
    Y_est = nan(size(k,1), n_obs*ny);

    % Find and count number of MKF filters
    n_obs_mkf = 0;
    n_obs_mkf_sp = 0;
    obs_mkf = nan(1, n_obs_mkf);
    obs_mkf_sp = nan(1, n_obs_mkf_sp);
    for i = 1:n_obs
        if ismember(observers{i}.type, ["MKF", "MKF_BM", "MKF_SF", ...
                "MKF_SF_DI", "MKF_SF_RODD", "MKF_SF_RODD95", ...
                "MKF_SP", "MKF_SF_DI", "MKF_SP_RODD"])
            n_obs_mkf = n_obs_mkf + 1;
            obs_mkf(n_obs_mkf) = i;
        end
        if ismember(observers{i}.type, ["MKF_SP", "MKF_SP_DI", ...
                "MKF_SP_RODD"])
            n_obs_mkf_sp = n_obs_mkf_sp + 1;
            obs_mkf_sp(n_obs_mkf_sp) = i;
        end
    end

    % Cell arrays to store data on each MKF observer
    MKF_i = cell(1, n_obs_mkf);
    MKF_X_est = cell(1, n_obs_mkf);
    MKF_trP = cell(1, n_obs_mkf);
    MKF_p_seq_g_Yk = cell(1, n_obs_mkf);
    for f = 1:n_obs_mkf
        obs = observers{obs_mkf(f)};
        MKF_i{f} = zeros(size(k));
        MKF_X_est{f} = nan(length(k), obs.nh*n);
        MKF_trP{f} = nan(length(k), obs.nh);
        MKF_p_seq_g_Yk{f} = nan(length(k), obs.nh);
    end
    MKF_SP_f_main = cell(1, n_obs_mkf_sp);
    MKF_SP_f_hold = cell(1, n_obs_mkf_sp);
    for f = 1:n_obs_mkf_sp
        obs = observers{obs_mkf_sp(f)};
        MKF_SP_f_main{f} = int16(zeros(length(k), obs.n_main));
        MKF_SP_f_hold{f} = int16(zeros(length(k), obs.n_hold));
    end
    % Initialize process
    xk = x0;

    % Choose debugging point
    k_stop = -1;  % no debugging
    %k_stop = 452;
    show_plots = false;

    % Start simulation at k = 0 (i = 1)
    for ki = 0:nT

        i = ki + 1;  % array index starts at 1

        % Use this for debugging - set k_stop above
        if ki == k_stop
            show_plots = true;
            k_stop = ki + 1;
        else
            show_plots = false;
        end

        % Get current input
        uk = U(i, :)';

        % Calculate current process outputs
        yk = meas_fcn(xk, uk, Ts, params);
        yk_m = yk + V(i,:)';

        % Record system states and outputs
        X(i, :) = xk';
        Y(i, :) = yk';
        Y_m(i, :) = yk_m';

        % Record observer estimates from previous timestep
%         for j = 1:n_obs
%             Xkp1_est(i, n*(j-1)+1:n*j) = observers{j}.xkp1_est';
%             Pkp1(i, n*(j-1)+1:n*j) = diag(observers{j}.Pkp1);
%             Ykp1_est(i, ny*(j-1)+1:ny*j) = observers{j}.ykp1_est';
%         end

        % Input vector for observers with measured inputs
        uk_m = uk(u_known);

        for j = 1:n_obs
            obs = observers{j};  % handle

            % Record observer predictions (this is over-written
            % below if the observer produces an updated estimate).
            if isprop(obs, 'Pkp1') || isfield(obs, 'Pkp1')
                X_est(i, n*(j-1)+1:n*j) = obs.xkp1_est';
                P(i, n*(j-1)+1:n*j) = diag(obs.Pkp1);
                if isprop(obs, 'ykp1_est') || isfield(obs, 'ykp1_est')
                    Y_est(i, ny*(j-1)+1:ny*j) = obs.ykp1_est';
                end
            end

            % Update observers using measurements uk_m and yk_m
            switch obs.type
                case {"EKF", "MEKF"}
                    % Update observer estimates
                    obs = update_MEKF(obs, yk_m, uk_m, Ts);
                    observers{j} = obs;  % copy updates to original
                otherwise
                    % Update observer estimates
                    obs.update(yk_m, uk_m);
            end
        end

        % Loop over each observer
        for j = 1:n_obs
            obs = observers{j};  % handle

            % Record observer updated estimates
            if isprop(obs, 'Pk') || isfield(obs, 'Pk')
                X_est(i, n*(j-1)+1:n*j) = obs.xk_est';
                P(i, n*(j-1)+1:n*j) = diag(obs.Pk);
                Y_est(i, ny*(j-1)+1:ny*j) = obs.yk_est';
            end

            % Record additional data for certain observers
            switch obs.type

                case {"MKF_S"}
                    % Multi-model Kalman filters and including MKF
                    % with sequence fusion

                    % Save simulation data for plotting later
                    f_mkf = find(obs_mkf == j);
                    MKF_i{f_mkf}(i) = obs.i;
                    for f = 1:obs.nh
                        MKF_X_est{f_mkf}(i, n*(f-1)+1:n*f) = ...
                            obs.filters.Xk_est(:,:,f)';
                        MKF_trP{f_mkf}(i, f) = trace(obs.filters.Pk(:,:,f));
                    end
                    MKF_p_seq_g_Yk{f_mkf}(i, 1:obs.nh) = obs.p_seq_g_Yk(1:obs.nh)';

                case {"MKF_BM", "MKF_GPB1", "MKF_GPB2", "MKF_SF", ...
                      "MKF_SF_DI", "MKF_SF_RODD", "MKF_SF_RODD95"}
                    % Multi-model Kalman filters with branching and
                    % merging. In this case we save the merged hypotheses
                    % not those associated with the filters.

                    % Save simulation data for plotting later
                    f_mkf = find(obs_mkf == j);
                    MKF_i{f_mkf}(i) = obs.i;
                    for f = 1:obs.nm
                        MKF_X_est{f_mkf}(i, n*(f-1)+1:n*f) = ...
                            obs.merged.Xk_est(:,:,f)';
                        MKF_trP{f_mkf}(i, f) = trace(obs.merged.Pk(:,:,f));
                    end
                    MKF_p_seq_g_Yk{f_mkf}(i, 1:obs.nm) = obs.merged.p_seq_g_Yk';

                case {"MKF_SP", "MKF_SP_DI", "MKF_SP_RODD"}
                    % MKF observer with sequence pruning

                    % Save simulation data for plotting later
                    % Save simulation data for plotting later
                    f_mkf = find(obs_mkf == j);
                    for f = 1:obs.nh
                        MKF_X_est{f_mkf}(i, n*(f-1)+1:n*f) = ...
                            obs.filters.Xkp1_est(:,:,f)';
                        MKF_trP{f_mkf}(i, f) = trace(obs.filters.Pkp1(:,:,f));
                    end
                    MKF_p_seq_g_Yk{f_mkf}(i, :) = obs.p_seq_g_Yk';

                    f_mkf_sp = find(obs_mkf_sp == j);
                    MKF_SP_f_main{f_mkf_sp}(i, :) = obs.f_main;
                    MKF_SP_f_hold{f_mkf_sp}(i, :) = obs.f_hold;

                case "EKF"
                    % Extended Kalman filters
                    error("NotImplementedError")

                case "MEKF"
                    % Extended multi-model Kalman filters
                    error("NotImplementedError")
            end
        end

        % Use this for step-by-step debugging
        if show_plots
            fprintf("%d: t = %g\n", ki, t(i));
            disp("k_stop reached");  % debug here
            show_plots = false;
        end

        % Calculate process states in next timestep
        wk = W(i, :)';  % process noise
        xk = state_fcn(xk, uk, Ts, params) + wk;

    end

    % Record final observer estimates
    if isprop(obs, 'Pk') || isfield(obs, 'Pk')
        X_est(i, n*(j-1)+1:n*j) = obs.xk_est';
        P(i, n*(j-1)+1:n*j) = diag(obs.Pk);
        Y_est(i, ny*(j-1)+1:ny*j) = obs.yk_est';
    end

    % Calculate output estimation errors
    E_obs = repmat(Y, 1, n_obs) - Y_est;

    % Main simulation results table
    sim_out.data = [table(k,t) ...
        table(U,V,W,X,X_est,P,Y,Y_m,Y_est,E_obs) ...
        input_data(:, other_data_labels)];

    % Data on MKF filters
    sim_out.MKF_obs = obs_mkf;
    sim_out.MKF_i = MKF_i;
    sim_out.MKF_X_est = MKF_X_est;
    sim_out.MKF_trP = MKF_trP;
    sim_out.MKF_p_seq_g_Yk = MKF_p_seq_g_Yk;
    sim_out.MKF_SP_f_main = MKF_SP_f_main;
    sim_out.MKF_SP_f_hold = MKF_SP_f_hold;
    sim_out.observers = observers;

end