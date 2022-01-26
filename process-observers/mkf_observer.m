function obs = mkf_observer(A,B,C,D,Ts,P0,Q,R,seq,T,d,label,x0)
% obs = mkf_observer(A,B,C,D,Ts,P0,Q,R,seq,T,d,label,x0)
%
% Creates a struct for simulating a multi-model Kalman filter
% for state estimation of a Markov jump linear system.
%
% Arguments:
%	A, B, C, D : cell arrays containing discrete-time system
%       matrices for each switching system modelled.
%   Ts : sample period.
%   P0 : cell array of initial covariance matrices of the 
%       state estimates for each filter.
%   Q : cell array of process noise covariance matrices for
%       each switching system.
%   R : cell array of output measurement noise covariance
%       matrices for each switching system.
%   seq : model indicator sequences for each filter (in rows).
%   T : transition probabity matrix of the Markov switching
%       process.
%   d : detection interval length in number of sample periods.
%   label : string name.
%   x0 : intial state estimates (optional).
%

    % Number of switching systems
    nj = numel(A);

    % System dimensions
    [n, nu, ny] = check_dimensions(A{1}, B{1}, C{1}, D{1});
    if nargin < 13
        x0 = zeros(n,1);
    end
    obs.A = A;
    obs.B = B;
    obs.C = C;
    obs.D = D;
    obs.Ts = Ts;
    obs.Q = Q;
    obs.R = R;
    obs.seq = seq;
    obs.T = T;
    obs.d = d;
    obs.label = label;

    % Check transition probability matrix
    assert(all(abs(sum(T, 2) - 1) < 1e-15), "ValueError: T")

    % Initialize covariance matrix of estimation errors
    obs.P = P0;

    % Check all other system matrix dimensions have same 
    % input/output dimensions and number of states.
    for j = 2:nj
        [n_j, nu_j, ny_j] = check_dimensions(A{j}, B{j}, C{j}, D{j});
        assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
            "ValueError: A, B, C, D")
    end

    % Number of filters required
    n_filt = size(seq, 1);

    % Fusion horizon length
    f = size(cell2mat(seq), 2);

    % Sequence index and counter for prob. updates
    % obs.i(1) is the sequence index (1 <= i(1) <= obs.f)
    % obs.i(2) is the counter for prob. updates (1 <= i(2) <= obs.d)
    obs.i = nan(1, 2);
    obs.i_next = int16([1 1]);

    % Initialize conditional probabilities: all equally likely
    obs.p_seq_g_Yk = ones(n_filt, 1) ./ n_filt;

    % Empty vectors to store values for filter calculations
    % p(y(k)|Gamma(k),Y(k-1))
    obs.p_yk_g_seq_Ykm1 = zeros(n_filt, 1);
    % Pr(gamma(k)|Y(k-1))
    obs.p_gammak_g_Ykm1 = zeros(n_filt, 1);
    % Pr(gamma(k))
    obs.p_gamma_k = zeros(n_filt, 1);
    % Pr(Gamma(k)|Y(k-1))
    obs.p_seq_g_Ykm1 = zeros(n_filt, 1);

    % Create multi-model observer
    obs.filters = cell(n_filt, 1);
    fmt = strcat('%s%0', char(string(strlength(string(n_filt)))), 'd');
    for i = 1:n_filt
        label_i = sprintf(fmt,label,i);
        % Initialize each filter with system #1
        obs.filters{i} = kalman_filter(A{1},B{1},C{1},D{1},Ts,P0{i}, ...
            Q{1},R{1},label_i,x0);
    end

    % Initialize estimates
    obs.xkp1_est = x0;
    obs.ykp1_est = obs.C{1} * obs.xkp1_est;

    % Save useful variables in struct
    obs.nj = nj;
    obs.n = n;
    obs.nu = nu;
    obs.ny = ny;
    obs.n_filt = n_filt;
    obs.f = f;

end