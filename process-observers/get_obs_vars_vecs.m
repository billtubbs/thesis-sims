function varargout = get_obs_vars_vecs(obs)
% varargout = get_obs_vars_vecs(obs)
% Returns data vectors containing values for all the
% time-varying variables of the observer. This is used
% in S-functions for storing the variable data in the
% Dwork memory objects allocated by Simulink.
%
% The number of data vectors returned depends on the
% type of observer. Some observers have integer
% variables which need to be stored in a separate
% Dwork vector.
%
% See function set_obs_vars_vecs for the reverse 
% operation - i.e. creating the value vectors from the
% observer struct.
%
% Examples
% >> vec_double = get_obs_vars_vecs(obs1);
% >> [vec_double, vec_int16] = get_obs_vars_vecs(obs2);
%

% Get observer variables
vars = get_obs_vars(obs);

switch obs.type

    case {"KFPSS", "LB"}  % steady-state filters

        vars_double = {vars.xkp1_est, vars.ykp1_est};

        % Convert variables to vectors
        vdata = make_data_vectors(vars_double);
        varargout{1} = cell2mat(vdata.vecs);

    case "KFFSS"  % steady-state Kalman filter - filtering form

        vars_double = {vars.xkp1_est};

        % Convert variables to vectors
        vdata = make_data_vectors(vars_double);
        varargout{1} = cell2mat(vdata.vecs);

    case "KFP"  % Kalman filter - prediction form

        vars_double = {vars.xkp1_est, vars.ykp1_est, vars.Pkp1};

        % Convert variables to vectors
        vdata = make_data_vectors(vars_double);
        varargout{1} = cell2mat(vdata.vecs);

    case {"KFF", "SKF"}  % Kalman filters and switching KF

        vars_double = {vars.xkp1_est, vars.Pkp1};

        % Convert variables to vectors
        vdata = make_data_vectors(vars_double);
        varargout{1} = cell2mat(vdata.vecs);

    case {"SKF_S"}  % Kalman filters and switching KF

        vars_double = {vars.xkp1_est, vars.Pkp1};

        % Convert variables to vectors
        vdata = make_data_vectors(vars_double);
        vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');
        varargout{1} = cell2mat(vdata.vecs);
        varargout{2} = cell2mat(vdata_int16.vecs);

    case {"MKF", "MKF_S", "MKF_SF", "MKF_SF_RODD95"}  % multi-model observers

        vars_double = {vars.xkp1_est, vars.p_seq_g_Yk, ...
            vars.xkp1_est_f, vars.Pkp1_f};

        % Convert variables to vectors
        vdata = make_data_vectors(vars_double);
        vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');
        varargout{1} = cell2mat(vdata.vecs);
        varargout{2} = cell2mat(vdata_int16.vecs);

    case "MKF_SF_RODD"

        vars_double = {vars.xkp1_est, vars.p_seq_g_Yk, ...
            vars.xkp1_est_f, vars.Pkp1_f};

        % Convert variables to vectors
        vdata = make_data_vectors(vars_double);
        vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');
        varargout{1} = cell2mat(vdata.vecs);
        varargout{2} = cell2mat(vdata_int16.vecs);

    case {"MKF_SP", "MKF_SP_DI", "MKF_SP_RODD"}

        vars_double = {vars.xkp1_est, vars.p_seq_g_Yk, ...
            vars.xkp1_est_f, vars.Pkp1_f};

        % Convert variables to vectors
        vdata = make_data_vectors(vars_double);
        vdata_int16 = make_data_vectors(struct2cell(vars.int16)', 'int16');
        varargout{1} = cell2mat(vdata.vecs);
        varargout{2} = cell2mat(vdata_int16.vecs);

    otherwise
        error('Value error: observer type not recognized')

end