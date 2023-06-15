function [idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, ...
    nj, nf)
% [idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj, nf)
% Generates three sequences of vectors of indices used in
% the sequence fusion sub-optimal multi-model observer 
% algorithm described by Robertson and Lee (1998). These 
% vectors determine the branching, mode transitions, and 
% sequence merging steps for a given set of sequences defined 
% over the fusion horizon. This version of the calculations
% differs from seq_fusion_indices_inv(seq, nj) in that it
% finds a set of branching, transition and merging steps
% for each time-step over the horizon, which may not be
% the same.
%
% For more details, see eqn's 22 and 23 in the Robertson and 
% Lee paper referenced below.
%
% Arguments:
%   seq : (nh, ns) integer double
%       Set of nh mode sequences to be modelled over the
%       fustion horizon of length f. At every time
%       instant, past sequences are merged into one of
%       these sequences and others are discarded.
%   nj : integer double
%       Number of system modes (2 for a single RODD).
%   nf : integer double (optional)
%       Fusion horizon. If not provided, nf is set to
%       the length of the sequence (i.e. size(seq, 2)).
%
% Returns:
%   idx_branch : (nb, 1) integer double
%       This column vector determines how the nh hypotheses
%       will be branched (i.e. split) into nb new hypotheses
%       at the next time instant.
%   idx_modes : (nb, 1) integer double
%       This column vector indicates the mode transitions to
%       extend the nb new hypotheses to the next time instant.
%   idx_merge : (nb, 1) integer double
%       This column vector determines how the nb hypotheses
%       will be merged (i.e. combined) into nh hypotheses
%       at the next time instant.
%
% Example 1:
% >> nj = 2;  % Number of system modes
% >> f = 3;  % Fusion horizon
% >> m = 1;  % Maximum number of shocks
% >> seq = cell2mat(shock_combinations_lte(f, m))
% 
% seq =
% 
%   4×3 int16 matrix
% 
%    1   1   1
%    2   1   1
%    1   2   1
%    1   1   2
% 
% 
% >> [idx_branch, idx_modes, idx_merge] = seq_fusion_indices(seq, nj);
% >> cell2mat(idx_branch)
% 
% ans =
% 
%      1     1     1
%      1     1     1
%      2     2     2
%      2     3     3
%      3     3     4
%      4     4     4
% 
% >> cell2mat(idx_modes)
% 
% ans =
% 
%      1     1     1
%      2     2     2
%      1     1     1
%      2     1     1
%      1     2     1
%      1     1     2
% 
% >> cell2mat(idx_merge)
% 
% ans =
% 
%      1     1     1
%      2     3     4
%      1     2     2
%      2     1     3
%      3     3     1
%      4     4     4
% 
% 
% References:
%  -  Robertson, D. G., & Lee, J. H. (1998). A method for the
%     estimation of infrequent abrupt changes in nonlinear 
%     systems. Automatica, 34(2), 261–270. 
%     https://doi.org/10.1016/S0005-1098(97)00192-1
%
    arguments
        seq double
        nj {mustBeInteger}
        nf {mustBeInteger} = size(seq, 2)
    end

    % Number of merged hypotheses modelled
    nh = size(seq, 1);

    % Length of sequences
    ns = size(seq, 2);

    % Cell arrays to hold results
    idx_branch = cell(1, nf);
    idx_modes = cell(1, nf); 
    idx_merge = cell(1, nf);

    % Make copy of sequence starting in first time instant
    seq_i = seq;

    for i = 1:ns

        % Sequence starting at next time instant
        seq_ip1 = circshift(seq_i, -1, 2);

        % Construct indices for all possible branching and 
        % mode transition steps for next step of the sequence
        idx_branch{i} = reshape(repmat((1:nh), nj, 1), [], 1);
        idx_modes{i} = reshape(repmat((1:nj)', 1, nh), [], 1);

        % Hypothesis sequences from time k-f to k:
        seq_kmf_to_k = [seq_i(idx_branch{i}, end-nf+1:end) idx_modes{i}];

        % Drop first (oldest) values in sequences to get
        % hypothesis sequences from time k-f+1 to k:
        seq_kmfp1_to_k = seq_kmf_to_k(:, 2:end);

        % Drop sequences not defined to be modelled at next
        % time instant (i.e. not found in seq)
        seq_to_keep = ismember(seq_kmfp1_to_k, seq_ip1(:, end-nf+1:end), 'rows');
        seq_kmfp1_to_k = seq_kmfp1_to_k(seq_to_keep, :);
        idx_branch{i} = idx_branch{i}(seq_to_keep, :);
        idx_modes{i} = idx_modes{i}(seq_to_keep, :);

        % Check re-merged sequences match the next sequence
        assert(isequal(unique(seq_kmfp1_to_k,'rows'), ...
            sortrows(seq_ip1(:, end-nf+1:end))))

        % Index of identical sequences (rows of seq_kmfp1_to_k)
        % which should be merged
        [~, idx_merge{i}] = ismember(seq_kmfp1_to_k, ...
            seq_ip1(:, end-nf+1:end), 'rows');

        seq_i = seq_ip1;

    end

end