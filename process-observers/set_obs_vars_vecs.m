function obs = set_obs_vars_vecs(obs, varargin)
% obs = set_obs_vars_vecs(obs, varargin)
% Takes a certain number of vectors containing values
% for variables and sets those variables in the observer 
% struct. This is used in S-functions for 'unpacking'
% the variable data from the Dwork memory objects
% allocated by Simulink.
%
% See function get_obs_vars_vecs for the reverse 
% operation - i.e. creating the value vectors from the 
% observer struct.
%
% Example
% >> obs = set_obs_vars_vecs(obs, vec_double, vec_int16);
%
% In this example, vec_double contains the appropriate number
% of real double values and vec_int16 contains the
% appropriate number of integer values to specify all the
% variables in obs.
%

    n = obs.n;
    ny = obs.ny;

    switch obs.type

        case {"LB", "KFPSS"}  % steady-state filters

            assert(nargin == 2)
            vec_double = varargin{1};

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. ykp1_est : size(ny, 1)

            vdata.types = {'double', 'double'};
            vdata.dims = {[n 1], [ny 1]};
            vdata.n_els = {n, ny};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));

            % Unpack data vectors - doubles
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.ykp1_est = vars_double{2};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case "KFFSS"  % steady-state Kalman filter - filtering form

            assert(nargin == 2)
            vec_double = varargin{1};

            % double vector contents:
            % 1. xkp1_est : size(n, 1)

            vdata.types = {'double'};
            vdata.dims = {[n 1]};
            vdata.n_els = {n};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));

            % Unpack data vectors - doubles
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case "KFP"  % Kalman filter - prediction form

            assert(nargin == 2)
            vec_double = varargin{1};

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. ykp1_est : size(ny, 1)
            % 3. Pkp1 : size(n, n)

            vdata.types = {'double', 'double', 'double'};
            vdata.dims = {[n 1], [ny 1], [n n]};
            vdata.n_els = {n, ny, n*n};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));

            % Unpack data vectors - doubles
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.ykp1_est = vars_double{2};
            vars.Pkp1 = vars_double{3};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case {"KFF", "SKF"}  % Kalman filters and switching KF

            assert(nargin == 2)
            vec_double = varargin{1};

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. Pkp1 : size(n, n)

            vdata.types = {'double', 'double'};
            vdata.dims = {[n 1], [n n]};
            vdata.n_els = {n, n*n};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));

            % Unpack data vectors - doubles
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.Pkp1 = vars_double{2};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

%         case "SKF"  % Switching Kalman filter
% 
%             assert(nargin == 3)
%             vec_double = varargin{1};
%             vec_int16 = varargin{2};
% 
%             % double vector contents:
%             % 1. xkp1_est : size(n, 1)
%             % 2. Pkp1 : size(n, n)
%             % 3. rk : size(1, 1)
% 
%             % int16 vector contents:
%             % 1. rk : size(1, 1)
% 
%             vdata.types = {'double', 'double', 'double'};
%             vdata.dims = {[n 1], [n n], [1 1]};
%             vdata.n_els = {n, n*n, 1};
%             vdata_int16.types = {'int16'};
%             vdata_int16.dims = {[1 1]};
%             vdata_int16.n_els = {1};
% 
%             % Add variables data
%             vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
%             vdata_int16.vecs = mat2cell(vec_int16', 1, ...
%                 cell2mat(vdata_int16.n_els));
% 
%             % Unpack data vectors - doubles
%             vars_double = unpack_data_vectors(vdata);
%             vars = struct();
%             vars.xkp1_est = vars_double{1};
%             vars.Pk = vars_double{2};
%             vars.rk = vars_double{3};
% 
%             % Unpack data vectors - integers
%             vars_int16 = unpack_data_vectors(vdata_int16);
%             vars.int16.rk = vars_int16{1};
% 
%             % Set all dynamic variables
%             obs = set_obs_vars(obs, vars);

        case "SKF_S"  % Switching Kalman filter with schedule

            assert(nargin == 3)
            vec_double = varargin{1};
            vec_int16 = varargin{2};

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. Pkp1 : size(n, n)

            % int16 vector contents:
            % 1. i : size(1, 1)
            % 2. i_next : size(1, 1)

            % Static data to unpack vectors
            vdata.types = {'double', 'double'};
            vdata.dims = {[n 1], [n n]};
            vdata.n_els = {n, n*n};
            vdata_int16.types = {'int16', 'int16'};
            vdata_int16.dims = {[1 1], [1 1]};
            vdata_int16.n_els = {1, 1};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
            vdata_int16.vecs = mat2cell(vec_int16', 1, ...
                cell2mat(vdata_int16.n_els));

            % Unpack data vectors - doubles
            % TODO: Maybe this is unnecessary, just set the obs
            % attributes directly.
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.Pkp1 = vars_double{2};

            % Unpack data vectors - integers
            vars_int16 = unpack_data_vectors(vdata_int16);
            vars.int16.i = vars_int16{1};
            vars.int16.i_next = vars_int16{2};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case "MKF"

            assert(nargin == 3)
            vec_double = varargin{1};
            vec_int16 = varargin{2};

            nh = obs.nh;

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. p_seq_g_Yk : size(nh, 1)
            % 3. xkp1_est_f (xkp1_est for each KF) : cell(1, nh)
            % 4. Pkp1_f (Pkp1 for each KF) : cell(1, nh)

            % int16 vector contents:
            % 1. rk : size(nh, 1)

            % Static data to unpack vectors
            vdata.types = {'double', 'double', 'double', 'double'};
            vdata.dims = {[n 1], [nh 1], [n 1 nh], [n n nh]};
            vdata.n_els = {n, nh, nh*n, nh*n*n};
            vdata_int16.types = {'int16'};
            vdata_int16.dims = {[nh 1]};
            vdata_int16.n_els = {nh};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
            vdata_int16.vecs = mat2cell(vec_int16', 1, ...
                cell2mat(vdata_int16.n_els));

            % Unpack data vectors - doubles
            % TODO: Maybe this is unnecessary, just set the obs
            % attributes directly.
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.p_seq_g_Yk = vars_double{2};
            vars.xkp1_est_f = vars_double{3};
            vars.Pkp1_f = vars_double{4};

            % Unpack data vectors - integers
            vars_int16 = unpack_data_vectors(vdata_int16);
            vars.int16.rk = vars_int16{1};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case {"MKF_S", "MKF_SF", "MKF_SF_RODD95"}

            assert(nargin == 3)
            vec_double = varargin{1};
            vec_int16 = varargin{2};

            nh = obs.nh;

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. p_seq_g_Yk : size(nh, 1)
            % 3. xkp1_est_f (xkp1_est for each KF) : cell(1, nh)
            % 4. Pkp1_f (Pkp1 for each KF) : cell(1, nh)

            % int16 vector contents:
            % 1. i : size(1, 1)
            % 2. i_next : size(1, 1)

            % Static data to unpack vectors
            vdata.types = {'double', 'double', 'double', 'double'};
            vdata.dims = {[n 1], [nh 1], [n 1 nh], [n n nh]};
            vdata.n_els = {n, nh, nh*n, nh*n*n};
            vdata_int16.types = {'int16', 'int16'};
            vdata_int16.dims = {[1 1], [1 1]};
            vdata_int16.n_els = {1, 1};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
            vdata_int16.vecs = mat2cell(vec_int16', 1, ...
                cell2mat(vdata_int16.n_els));

            % Unpack data vectors - doubles
            % TODO: Maybe this is unnecessary, just set the obs
            % attributes directly.
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.p_seq_g_Yk = vars_double{2};
            vars.xkp1_est_f = vars_double{3};
            vars.Pkp1_f = vars_double{4};

            % Unpack data vectors - integers
            vars_int16 = unpack_data_vectors(vdata_int16);
            vars.int16.i = vars_int16{1};
            vars.int16.i_next = vars_int16{2};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case {"MKF_SF_RODD"}

            assert(nargin == 3)
            vec_double = varargin{1};
            vec_int16 = varargin{2};

            nh = obs.nh;

            % {vars.xkp1_est, vars.p_seq_g_Yk, ...
            %  vars.rk, vars.xkp1_est_f, vars.Pkp1_f}

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. p_seq_g_Yk : size(nh, 1)
            % 3. xkp1_est_f (xkp1_est for each KF) : cell(1, nh)
            % 4. Pkp1_f (Pkp1 for each KF) : cell(1, nh)

            % int16 vector contents:
            % 1. rk : size(nh, 1)
            % 2. i : size(1, 1)
            % 3. i_next : size(1, 1)
            % 4. id : size(1, 1)
            % 5. id_next : size(1, 1)

            % Static data to unpack vectors
            vdata.types = {'double', 'double', 'double', 'double'};
            vdata.dims = {[n 1], [nh 1], [n 1 nh], [n n nh]};
            vdata.n_els = {n, nh, nh*n, nh*n*n};
            vdata_int16.types = {'int16', 'int16', 'int16', 'int16', 'int16'};
            vdata_int16.dims = {[nh 1], [1 1], [1 1], [1 1], [1 1]};
            vdata_int16.n_els = {nh, 1, 1, 1, 1};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
            vdata_int16.vecs = mat2cell(vec_int16', 1, ...
                cell2mat(vdata_int16.n_els));

            % Unpack data vectors - doubles
            % TODO: Maybe this is unnecessary, just set the obs
            % attributes directly.
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.p_seq_g_Yk = vars_double{2};
            vars.xkp1_est_f = vars_double{3};
            vars.Pkp1_f = vars_double{4};

            % Unpack data vectors - integers
            vars_int16 = unpack_data_vectors(vdata_int16);
            vars.int16.rk = vars_int16{1};
            vars.int16.i = vars_int16{2};
            vars.int16.i_next = vars_int16{3};
            vars.int16.id = vars_int16{4};
            vars.int16.id_next = vars_int16{5};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case {"MKF_SP"}

            assert(nargin == 3)
            vec_double = varargin{1};
            vec_int16 = varargin{2};

            nh = obs.nh;

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. p_seq_g_Yk : size(nh, 1)
            % 3. xkp1_est for each KF : cell(1, nh)
            % 4. Pkp1 for each KF : cell(1, nh)

            % int16 vector contents:
            % 1. rk : size(nh, 1)
            % 2. f_main : size(1, n_main)
            % 3. f_hold : size(1, n_hold)

            % Static data to unpack vectors
            vdata.types = {'double', 'double', 'double', 'double'};
            vdata.dims = {[n 1], [nh 1], [n 1 nh], [n n nh]};
            vdata.n_els = {n, nh, nh*n, nh*n*n};
            vdata_int16.types = {'int16', 'int16', 'int16'};
            vdata_int16.dims = {[nh 1], [1 obs.n_main], [1 obs.n_hold]};
            % Note all elements of this cell array must be of same class
            vdata_int16.n_els = {nh, obs.n_main, obs.n_hold};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
            vdata_int16.vecs = mat2cell(vec_int16', 1, ...
                cell2mat(vdata_int16.n_els));

            % Unpack data vectors - doubles
            % TODO: Maybe this is unnecessary, just set the obs
            % attributes directly.
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.p_seq_g_Yk = vars_double{2};
            vars.xkp1_est_f = vars_double{3};
            vars.Pkp1_f = vars_double{4};

            % Unpack data vectors - integers
            vars_int16 = unpack_data_vectors(vdata_int16);
            vars.int16.rk = vars_int16{1};
            vars.int16.f_main = vars_int16{2};
            vars.int16.f_hold = vars_int16{3};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        case {"MKF_SP_DI", "MKF_SP_RODD"}

            assert(nargin == 3)
            vec_double = varargin{1};
            vec_int16 = varargin{2};

            nh = obs.nh;

            % double vector contents:
            % 1. xkp1_est : size(n, 1)
            % 2. p_seq_g_Yk : size(nh, 1)
            % 3. xkp1_est for each KF : cell(1, nh)
            % 4. Pkp1 for each KF : cell(1, nh)

            % int16 vector contents:
            % 1. rk : size(nh, 1)
            % 2. id : size(1, 1)
            % 3. id_next : size(1, 1)
            % 4. f_main : size(1, n_main)
            % 5. f_hold : size(1, n_hold)

            % Static data to unpack vectors
            vdata.types = {'double', 'double', 'double', 'double'};
            vdata.dims = {[n 1], [nh 1], [n 1 nh], [n n nh]};
            vdata.n_els = {n, nh, nh*n, nh*n*n};
            vdata_int16.types = {'int16', 'int16', 'int16', 'int16', ...
                'int16'};
            vdata_int16.dims = {[nh 1], [1 1], [1 1], ...
                [1 obs.n_main], [1 obs.n_hold]};
            % Note all elements of this cell array must be of same class
            vdata_int16.n_els = {nh, 1, 1, obs.n_main, obs.n_hold};

            % Add variables data
            vdata.vecs = mat2cell(vec_double', 1, cell2mat(vdata.n_els));
            vdata_int16.vecs = mat2cell(vec_int16', 1, ...
                cell2mat(vdata_int16.n_els));

            % Unpack data vectors - doubles
            % TODO: Maybe this is unnecessary, just set the obs
            % attributes directly.
            vars_double = unpack_data_vectors(vdata);
            vars = struct();
            vars.xkp1_est = vars_double{1};
            vars.p_seq_g_Yk = vars_double{2};
            vars.xkp1_est_f = vars_double{3};
            vars.Pkp1_f = vars_double{4};

            % Unpack data vectors - integers
            vars_int16 = unpack_data_vectors(vdata_int16);
            vars.int16.rk = vars_int16{1};
            vars.int16.id = vars_int16{2};
            vars.int16.id_next= vars_int16{3};
            vars.int16.f_main = vars_int16{4};
            vars.int16.f_hold = vars_int16{5};

            % Set all dynamic variables
            obs = set_obs_vars(obs, vars);

        otherwise
            error("Value error: observer type not recognized")

    end