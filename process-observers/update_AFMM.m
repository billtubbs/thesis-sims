function obs = update_AFMM(obs, uk, yk)
% obs = update_AFMM(obs, uk, yk) updates the multi-
% model Kalman filter using the adaptive forgetting 
% through multiple models (AFMM) algorithm and 
% calculates the estimates of the states and output 
% at the next sample time.
%
% Sequence pruning method as described in Eriksson and
% Isaksson (1996):
% - At the updates, let only the most probable sequence,
%   i.e. with the largest weight of all the sequences
%   split into 2 branches.
% - After the update is completed, remove the sequence
%   with the smallest weight.
% - Renormalize the remaining weights to have unit sum.
%
% Restriction to above rules:
% - Do not cut branches immediately after they are born.
%   Let there be a certain minimum life length for all 
%   branches.
%
% Arguments:
%   obs : struct containing the multi-model Kalman filter
%       variables (see function mkf_filter).
%   uk : vector (nu, 1) of system inputs at the current 
%       sample time.
%   yk : vector (ny, 1) of system output measurements
%       at the current sample time.
%
% Reference:
%  -  Eriksson, P.-G., & Isaksson, A. J. (1996). Classification
%     of Infrequent Disturbances. IFAC Proceedings Volumes, 29(1), 
%     6614–6619. https://doi.org/10.1016/S1474-6670(17)58744-3
%

    % Implementation of filter pruning algorithm
    %
    % - As described in Eriksson and Isaksson (1996):
    %
    %   "Do not cut branches immediately after they are 
    %    born. Let there be a certain minimum life length 
    %    for all branches."
    %
    % To implement this rule, filters are organized into
    % two groups:
    % 
    %  1. Main group : obs.filters{f} for f in obs.f_hold
    %     Longer-surviving filters which can be killed off.
    %  2. Holding group : obs.filters{f} for f in obs.f_main
    %     Recently created filters which cannot be killed 
    %     off and are held for n_min time periods before
    %     being transferred to the main group.
    %
    % Note that for n_dist input disturbances, n_dist new
    % sequences are created each time step, so the
    % holding group needs to be n_dist*obs.n_min in size.

    % Index of current most likely sequence
    [~, f_max] = max(obs.p_seq_g_Yk);

    % Set next sequence value to 0 (no shock) for all
    % sequences
    for f = 1:obs.n_filt
        obs.seq{f}(:, obs.i_next(1)) = 0;
    end

    % Number of disturbances
    nw = size(obs.epsilon, 1);

    % Consistency checks - can be removed later
    assert(size(obs.f_hold, 2) == obs.n_hold)
    assert(size(obs.f_main, 2) == obs.n_main)
    assert(size(obs.f_unused, 2) == obs.n_filt-1)
    comb = [obs.f_hold obs.f_main obs.f_unused];
    comb = sort(nonzeros(comb))';
    assert(isequal(comb, 1:obs.n_filt))

    % Right-shift all filters in holding group. This causes
    % the last nw values to 'roll-over' to the left of f_hold.
    % e.g. circshift([1 2 3], 1) -> [3 1 2]
    obs.f_hold = circshift(obs.f_hold, nw);

    % Check if holding group is not yet full
    f_move = obs.f_hold(1:nw);
    i_move = f_move > 0;

    % Select filters to be moved out of holding group and
    % put into main group.
    if sum(i_move) > 0
        [obs.f_main, f_to_replace] = add_to_group_with_replacement(...
            obs.f_main, f_move(i_move), obs.p_seq_g_Yk);
    else
        f_to_replace = [];
    end

    % f_to_replace now contains all filter indices no longer in
    % either the main or holding group (i.e. those pruned due 
    % to their low likelihood). These are therefore now
    % available for use for new sequences.  However, during the
    % the first few time steps while f_main and f_hold are not
    % full, f_to_replace may be empty [].
    if numel(f_to_replace) < nw
        [~, obs.f_unused, f_to_replace] = ...
            take_from_2_groups(f_to_replace, obs.f_unused, nw);
    end
    obs.f_hold(1:nw) = f_to_replace;

    % Split most probable filter and add new filter(s)
    % to start of holding group
    for i = 1:nw
        label = obs.filters{f_to_replace(i)}.label;
        obs.filters{f_to_replace(i)} = obs.filters{f_max};
        obs.filters{f_to_replace(i)}.label = label;
        obs.p_seq_g_Yk(f_to_replace(i)) = obs.p_seq_g_Yk(f_max);
        obs.p_gamma_k(f_to_replace(i)) = obs.p_gamma_k(f_max);
        obs.p_seq_g_Ykm1(f_to_replace(i)) = obs.p_seq_g_Ykm1(f_max);
        obs.seq{f_to_replace(i)} = obs.seq{f_max};
        % Set next sequence value to 1 (shock)
        obs.seq{f_to_replace(i)}(:, obs.i_next(1)) = i;
    end

    % TODO: Add online noise variance estimation with
    % forgetting factor as described in Anderson (1995)
    % equation (4.3).

    % Run generic MKF updates and probability calcs first
    obs = update_MKF(obs, uk, yk);

end



