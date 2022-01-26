function obs = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon, ...
    sigma_wp,Q0,R,n_filt,f,n_min,label,x0)
% obs = mkf_observer_AFMM(A,B,C,D,Ts,u_meas,P0,epsilon,sigma_wp, ...
%     Q0,R,n_filt,f,n_min,label,x0)
%
% Creates a struct for simulating the multi-model observer
% using the adaptive forgetting through multiple models 
% (AFMM) method for state estimation in the presence of 
% infrequently-occurring deterministic disturbances, as 
% described in Eriksson and Isaksson (1996).
%
% Arguments:
%   A, B, C, D : matrices of the discrete time state-space
%       system representing the augmented system (including
%       disturbances and unmeasured inputs).
%   Ts : sample period.
%   u_meas : binary vector indicating measured inputs.
%   P0 : Initial value of covariance matrix of the state
%       estimates.
%   epsilon : probability of a shock disturbance.
%   sigma_wp : standard deviation of shock disturbances.
%   Q0 : Process noise covariance matrix (n, n) with 
%        variances for each state on the diagonal. The  
%        values for states impacted by the unmeasured
%        input disturbances should be set to zero as the
%        appropriate variances will be added by the
%        algorithm during observer updates.
%   R : output measurement noise covariance matrix (ny, ny).
%   n_filt : number of models (Kalman filters) to utilise.
%   f : length of disturbance sequences to record.
%   n_min : minimum life of cloned filters in number of
%       sample periods.
%   label : string name.
%   x0 : intial state estimates (optional).
%
% Reference:
%  -  Eriksson, P.-G., & Isaksson, A. J. (1996). Classification
%     of Infrequent Disturbances. IFAC Proceedings Volumes, 29(1), 
%     6614-6619. https://doi.org/10.1016/S1474-6670(17)58744-3
%

    % TODO: Expose spacing parameter d as an argument

    % Number of states
    n = check_dimensions(A, B, C, D);

    % Detection interval length in number of sample periods.
    d = 1;  % TODO: Make this a variable parameter

    % Initial state estimates
    if nargin == 15
        x0 = zeros(n,1);
    end

    % Observer model without disturbance noise input
    Bu = B(:, u_meas);
    Bw = B(:, ~u_meas);
    Du = D(:, u_meas);
    nw = sum(~u_meas);  % Number of input disturbances
    assert(nw > 0, "ValueError: u_meas");

    % Modified variances of shock signal over detection
    % interval (see (16) on p.264 of Robertson et al. 1998)
    var_wp = sigma_wp.^2 ./ d;

    % Construct process noise covariance matrices for each 
    % possible input disturbance
    Q = construct_Q(Q0, Bw, var_wp, u_meas);

    % Number of switching models
    nj = numel(Q);

    % Probabilities of no-shock, shock
    p_gamma = [1-epsilon epsilon]';

    if nw > 1

        % Possible combinations of each disturbance input:
        % Assume only one may occur in the same sample period
        Z = [zeros(1, nw); eye(nw)];
 
        % Modified indicator value probabilities
        p_gamma = prod(prob_gamma(Z', p_gamma), 1)';

        % Normalize so that sum(Pr(gamma(k))) = 1
        % TODO: Is this the right thing to do for sub-optimal approach?
        p_gamma = p_gamma ./ sum(p_gamma);

    end

    % Transition probability matrix
    % Note that for RODD disturbances Pr(gamma(k)) is
    % assumed to be an independent random variable.
    T = repmat(p_gamma', nj, 1);

    % Initialize indicator sequences
    seq = mat2cell(nan(n_filt, f), ones(1, n_filt), f);

    % Define filter groups ('main', 'holding' and 'unused')
    assert(n_min > 0)
    assert(n_filt > 0)
    n_hold = nw*n_min;
    n_main = n_filt - n_hold;
    f_hold = nan(1, n_hold);
    f_main = nan(1, n_main);
    f_main(1) = 1;  % initialize with one active filter
    f_unused = [2:n_filt];

    % Check there are enough filters in total to accommodate
    % those in the holding group + at least one in main group
    assert(n_main >= nw, "ValueError: n_filt is too low.")

    % System model doesn't change
    A = repmat({A}, 1, nj);
    Bu = repmat({Bu}, 1, nj);
    C = repmat({C}, 1, nj);
    Du = repmat({Du}, 1, nj);
    R = repmat({R}, 1, nj);

    % Initial covariance matrix is the same for all filters
    P0_init = repmat({P0}, 1, n_filt);

    % Create MKF observer struct
    obs = mkf_observer(A,Bu,C,Du,Ts,P0_init,Q,R,seq,T,d,label,x0);

    % AFMM is initialized differently
    obs.p_seq_g_Yk = [1; zeros(obs.n_filt-1, 1)];

    % Add additional variables used by AFMM observer
    obs.u_meas = u_meas;
    obs.f = f;
    obs.n_min = n_min;
    obs.n_main = n_main;
    obs.n_hold = n_hold;
    obs.f_main = f_main;
    obs.f_hold = f_hold;
    obs.f_unused = f_unused;
    obs.P0 = P0;
    obs.Q0 = Q0;
    obs.epsilon = epsilon;
    obs.sigma_wp = sigma_wp;
    obs.p_gamma = p_gamma;
    obs.nj = nj;

end


function Q = construct_Q(Q0, Bw, var_wp, u_meas)
% Q = construct_Q(Q0, B, var_wp, u_meas) 
% returns a cell array of different process noise 
% covariance matrices Qj for each filter model j = 
% 1:nj for the tracking of infrequently-occurring
% input disturbances.
%
% Arguments:
%   Q0 : nxn matrix containing variances for measured
%        states on the diagonal (all other elements 
%        are ignored).
%   Bw : system input matrix for the random shock signals
%       (n x nw).
%   var_wp : variances of shock disturbances over 
%       detection interval.
%   u_meas : binary vector indicating which inputs are
%        measured.
%

    % Number of states
    n = size(Bw, 1);

    % Check size of initial process covariance matrix
    assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

    % Number of input disturbances
    nw = sum(~u_meas);
    % Check size of disturbance input matrix
    assert(isequal(size(Bw), [n nw]))

    % TODO: This only works for disturbances with 2
    % states (i.e. shock/no shock). Could be extended
    % to other cases (e.g. no shock, small shock, big
    % shock)
    assert(size(var_wp, 2) == 2)

    % Number of switching models required. Assume only 
    % one shock per detection period is possible.
    nj = 1 + nw;

    % Array of variances of each shock combination
    var_x = repmat(diag(Q0), 1, nj);

    % Add noise variances for model 1 assuming no shocks
    % occurred
    var_x(:, 1) = var_x(:, 1) + Bw * var_wp(:, 1);

    % Add variances for models 2 to nj to reflect one
    % of each shocks occuring.
    idx = find(~u_meas);
    for i = 1:nw
        var_i = var_wp(:, 1);  % no shock
        var_i(i) = var_wp(i, 2);  % shock
        var_x(:, i+1) = var_x(:, i+1) + Bw * var_i;
    end

    Q = cell(1, nj);
    for j = 1:nj
        Q{j} = diag(var_x(:, j));
    end

end