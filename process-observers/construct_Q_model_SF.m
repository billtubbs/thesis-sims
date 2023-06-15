function [Q, p_rk, seq, p_seq] = construct_Q_model_SF(Q0, B, u_known, ...
    alpha, var_wp, n_di, m, nw)
% [Q, p_rk, seq, p_seq] = construct_Q_model_SF(Q0, B, u_known, ...
%     alpha, var_wp, n_di, m, nw)
% Constructs the parameters needed to model the sub-optimal
% multi-model algorithm using sequence fusion for the
% tracking of infrequently-occurring random disturbances.
% Returns a cell array of the process noise covariance 
% matrices Qj for each filter model j = 1:nj and a set of
% sequences of shock indicator values.
%
% Arguments:
%   Q0 : nxn matrix containing variances for measured
%        states on the diagonal (all other elements 
%        are ignored).
%   Bw : system input matrix for the random shock signals
%       (n x nw).
%   alpha : Probability of at least one shock in a 
%       detection interval.
%   var_wp : variances of shock disturbances over 
%       detection interval. See (16) on p.264 of Robertson 
%       et al. (1998).
%   n_di : number of detection intervals within the fusion 
%       horizon.
%   m : maximum number of disturbances over fusion horizon.
%   nw : number of independent input disturbances.
%

    % Generate indicator sequences
    seq = shock_combinations_lte(n_di*nw, m);

    % Input matrix for the unknown disturance inputs
    Bw = B(:, ~u_known);

    % Probabilities of no-shock / shock over detection interval
    % (this is named delta in Robertson et al. 1998)
    p_rk = [1-alpha'; alpha'];

    if nw == 1  % only one disturbance

        % Number of models (each with a different Q matrix)
        nj = 2;

        % Generate required Q matrices.
        Q = cell(1, nj);
        for i = 1:nj
            var_x = diag(Q0);
            var_x = var_x + Bw * var_wp(:, i);
            Q{i} = diag(var_x);
        end

    elseif nw > 1  % multiple disturbances

        % Note: In the case of more than one input disturbance
        % there may be multiple combinations of disturbances
        % occuring simultaneously. To simulate these, construct
        % a different Q matrix for each possible combination.

        % Reshape sequences into matrices with a row for each
        % input disturbance sequence
        seq = cellfun(@(x) reshape(x, nw, []), seq, ...
            'UniformOutput', false);

        % Find all unique combinations of simultaneous shocks
        [Z,~,ic] = unique(cell2mat(seq')', 'sorted', 'rows');

        % Number of models needed (each with a different Q matrix)
        nj = size(Z, 1);

        % Number of filters needed
        % TODO: arg f here is not actually the fusion horizon 
        % (which is f*d).  Should maybe use f/d here and assert no
        % remainder.
        n_filt = n_filters(m, n_di, nw);

        % Rearrange as one sequence for each filter and convert
        % back to cell array
        seq = reshape(ic', [], n_filt)';
        seq = mat2cell(seq, ones(n_filt, 1), n_di);

        % Generate required Q matrices
        Q = cell(1, nj);
        for i = 1:nj
            var_x = diag(Q0);
            ind = Z(i, :);
            % Modified variances of shock signal over detection
            % interval (see (16) on p.264 of Robertson et al. 1998)
            idx = sub2ind(size(var_wp), 1:nw, ind);
            var_x = var_x + Bw * var_wp(idx)';
            Q{i} = diag(var_x);
        end

        % Modify indicator value probabilities for
        % combinations of shocks
        p_rk = prod(prob_rk(Z', p_rk), 1)';

        % Normalize so that sum(Pr(gamma(k))) = 1
        % TODO: Is this the right thing to do for sub-optimal approach?
        p_rk = p_rk ./ sum(p_rk);

    end

    % Sequence probabilities Pr(Gamma(k))
    p_seq = prob_seq(seq, p_rk);

end