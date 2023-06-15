function vars = get_obs_vars(obs)
% vars = get_obs_vars(obs) returns a struct containing
% all the variables (i.e. time-varying properties) of the
% observer. The variables returned depends on the observer type.
%

% TODO: Write a test script for set_obs_vars and get_obs_vars
    switch(obs.type)

        case {"KFPSS", "LB"}  % Steady-state filters - prediction form

            % Get variables
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;

        case {"KFFSS"}  % Steady-state filter - filtering form

            % Get variables
            vars.xkp1_est = obs.xkp1_est;

        case {"KFP"}  % Kalman filter - prediction form

            % Get variables
            vars.xkp1_est = obs.xkp1_est;
            vars.ykp1_est = obs.ykp1_est;
            vars.Pkp1 = obs.Pkp1;

        case {"KFF", "SKF"}  % Kalman filters - filtering form

            % Get variables
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;

        case "SKF_S"  % Switching KF with schedule

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;

            % Get integer variables
            vars.int16.i = obs.i;
            vars.int16.i_next = obs.i_next;

        case "MKF"  % multi-model observers with sequences

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1;

            % Get integer variables
            vars.int16.rk = obs.rk;

        case {"MKF_S", "MKF_SF_RODD95"}  % multi-model observers with sequences

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1;

            % Get integer variables
            vars.int16.i = obs.i;
            vars.int16.i_next = obs.i_next;

        case {"MKF_SF_DI"}  % multi-model observers with detection intervals

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1;

            % Get integer variables
            vars.int16.rk = obs.rk;
            vars.int16.id = obs.id;
            vars.int16.id_next = obs.id_next;

        case {"MKF_SF_RODD"}  % multi-model observers with detection intervals

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1;

            % Get integer variables
            vars.int16.rk = obs.rk;
            vars.int16.i = obs.i;
            vars.int16.i_next = obs.i_next;
            vars.int16.id = obs.id;
            vars.int16.id_next = obs.id_next;

        case {"MKF_SFF", "MKF_SPF"}  % special multi-model Kalman filters

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.Pkp1 = obs.Pkp1;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = nan(obs.n, 1, obs.nh);
            vars.Pkp1_f = nan(obs.n, obs.n, obs.nh);
            for f = 1:obs.nh
               vars.xkp1_est_f(:,:,f) = obs.filters{f}.xkp1_est;
               vars.Pkp1_f(:,:,f) = obs.filters{f}.Pkp1;
            end

            % Get integer variables
            vars.int16.rk = obs.rk;
            vars.int16.id = obs.id;
            vars.int16.id_next = obs.id_next;
            if strcmp(obs.type, "MKF_SPF")
                % Additional variables used by adaptive sequence
                % pruning algorithms
                vars.int16.f_main = obs.f_main;
                vars.int16.f_hold = obs.f_hold;
                vars.int16.seq = obs.seq;
            end

        case {"MKF_SP"}  % multi-model observers

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1;

            % Get integer variables
            vars.int16.rk = obs.rk;
            vars.int16.f_main = obs.f_main;
            vars.int16.f_hold = obs.f_hold;

        case {"MKF_SP_DI", "MKF_SP_RODD"}  % multi-model observers

            % Get double variables
            vars.xkp1_est = obs.xkp1_est;
            vars.p_seq_g_Yk = obs.p_seq_g_Yk;
            vars.xkp1_est_f = obs.filters.Xkp1_est;
            vars.Pkp1_f = obs.filters.Pkp1;

            % Get integer variables
            vars.int16.rk = obs.rk;
            vars.int16.id = obs.id;
            vars.int16.id_next = obs.id_next;
            vars.int16.f_main = obs.f_main;
            vars.int16.f_hold = obs.f_hold;

        otherwise
            error("Value error: observer type not recognized")
    end

end