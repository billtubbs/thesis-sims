function [Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp)
% [Q, p_rk] = construct_Q_model_SP(Q0, B, u_known, alpha, sigma_wp);
% Constructs the parameters needed to model the sub-optimal
% multi-model algorithm using sequence pruning for the
% tracking of infrequently-occurring random disturbances.
% Returns a cell array of the process noise covariance
% matrices Qj for each filter model j = 1:nj for the 
% estimation of systems with infrequently-occurring input 
% disturbances.
%
% Arguments:
%   Q0 : (n, n)
%       Matrix containing variance values for process
%       states. Only values in the rows and columns corresponding 
%       to the process states are used, usually the upper left
%       block from (1, 1) to (n-nw, n-nw). The remaining
%       values corresponding to the covariances of the input 
%       disturbance model states are over-written at initialization.
%   B : (n, nu) double
%       System input matrix.
%   u_known : (nu, 1) logical
%       Boolean vector to indicate which inputs are known.
%   alpha : double
%       Probability of at least one shock in a detection interval.
%   sigma_wp : (1 nw) cell array
%       Standard deviations of disturbances. Each element of
%       the cell array is either a scalar for a standard (Gaussian)
%       noise or a (2 1) vector for a random shock disturbance.
%

    % Number of variances specified for each disturbance input
    ns = cellfun(@(s) size(s, 2), sigma_wp);

    % TODO: Currently, this function only works for disturbances
    % with up to 2 modes (e.g. shock, no shock). Could be extended
    % to other cases (e.g. no shock, small shock, big shock)
    assert(all(ns <= 2), "Only binary shocks supported")

    % Positions of unmeasured inputs in input signal
    idx_w = find(~u_known);

    % Number of unmeasured inputs (i.e. disturbances)
    nw = numel(idx_w);
    assert(numel(sigma_wp) == nw)

    % Number of switching noise models
    n_switch = sum(ns > 1);

    % Number of modes (= number of models required)
    % With this sequence pruning algorithm, no shocks may
    % occur simultaneously
    nj = 1 + sum(ns - 1);

    % Determine which disturbances have a switching variance
    % and construct all combinations of input variances.
    idx_switch = 1 + cumsum(ns-1);
    var_wp = zeros(nw, nj);
    for i = 1:nw
        if ns(i) == 1
            % Standard Gaussian noise input
            var_wp(i, :) = sigma_wp{i}^2;  % same for all modes
        else
            % Switching Gaussian noise input
            var_wp(i, :) = sigma_wp{i}(1)^2;
            var_wp(i, idx_switch(i)) = sigma_wp{i}(2)^2;  % different
        end
    end

    % Input matrix for the unknown disturance inputs
    Bw = B(:, ~u_known);

    % Generate the Q matrix for each system mode
    % Set all covariances on rows, cols of disturbance 
    % model to zero (i.i.d assumption).
    n = length(Q0);
    Q0(:, idx_w) = 0;
    Q0(idx_w, :) = 0;
    Q = cell(1, nj);
    for i = 1:nj
        Q{i} = Q0;
        % Scale variances using noise input matrix
        var_x = Bw * var_wp(:, i);
        % Add variances to elements of diagonal
        Q{i}(logical(eye(n))) = diag(Q{i}) + var_x;
    end

    if n_switch > 0
        % Probabilities of no-shock, shock
        p_rk = [1-alpha alpha]';

        if n_switch > 1

            % Possible combinations of each disturbance input:
            % Assume only one may occur in the same sample period
            Z = [ones(1, nw); eye(nw)+1];

            % Modified indicator value probabilities
            p_rk = prod(prob_rk(Z', p_rk), 1)';

            % Normalize so that sum(Pr(gamma(k))) = 1
            % TODO: Is this the right thing to do for sub-optimal approach?
            % No I don't think so. If hypotheses don't represent all
            % possible combinations then this is good to know. However,
            % the transition probability matrix must be normalized?
            % p_rk = p_rk ./ sum(p_rk);
        end
    else
        p_rk = 1;
    end
end