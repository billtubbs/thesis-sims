function [Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m, rk_m] = merge_estimates(Xk, Pk, Yk, ...
    p_seq_g_Yk, rk, idx)
% [Xk_m, Pk_m, p_seq_g_Yk_m] = merge_estimates(Xk, Pk, Yk, ...
%     p_seq_g_Yk, rk, idx)
% Merge (i.e. mix) a set of multiple-model Kalman filter
% estimates and probabilities according to the modes specified 
% in vector idx. Note that the sets of states, Xk, covariance 
% matrices, Pk, and output estimates, Yk, are 3-dimensional 
% arrays with the 3rd dimension represeting the hypotheses 
% over which to merge.
%
% Example:
% >> Xk = cat(3, 1, 2, 2, 3);
% >> Pk = cat(3, 10, 11, 12, 13);
% >> Yk = cat(3, 2, 4, 4, 6);
% >> p_seq_g_Yk = [0.8 0.2 0.2 0.8]';
% >> idx = [1 1 2 2]';
% >> [Xk_m,Pk_m,Yk_m,p_seq_g_Yk_m] = ...
%     merge_estimates(Xk,Pk,Yk,p_seq_g_Yk,rk,idx)
% 
% Xk_m(:,:,1) =
% 
%     1.2000
% 
% 
% Xk_m(:,:,2) =
% 
%     2.8000
% 
% 
% Pk_m(:,:,1) =
% 
%    10.3600
% 
% 
% Pk_m(:,:,2) =
% 
%    12.9600
% 
% 
% Yk_m(:,:,1) =
% 
%     2.4000
% 
% 
% Yk_m(:,:,2) =
% 
%     5.6000
% 
% 
% p_seq_g_Yk_m =
% 
%      1
%      1
%

    % Number of states
    n = size(Xk, 1);

    % Number of Outputs
    ny = size(Yk, 1);

    % Number of hypothesis before merge
    nh = size(p_seq_g_Yk, 1);

    % Add a small value to probabilities to prevent division
    % by zero. TODO: Is this a good idea?
    p_seq_g_Yk = p_seq_g_Yk + eps;
    p_seq_g_Yk = p_seq_g_Yk ./ sum(p_seq_g_Yk);

    % Weighting to apply to hypotheses
    weights = reshape(p_seq_g_Yk, 1, 1, nh);

    if nargin < 6
        % If rk and idx are not provided, do a complete merge
        % to one set of estimates

        % Merge state and output estimates
        Xk_m = sum(weights .* Xk, 3);
        if ~isempty(Yk)
            Yk_m = sum(weights .* Yk, 3);
        else
            Yk_m = [];
        end

        % Merge state estimation error covariance
        X_diffs = Xk_m - Xk;
        Pk_m = sum(weights .* (Pk + ...
            pagemtimes(X_diffs, pagetranspose(X_diffs))), 3);

        % Merged probabilities (sum to 1)
        p_seq_g_Yk_m = 1;

        % If rk_m is requested, return the mode with the highest
        % probability
        if nargout == 5
            rk_m = rk(p_seq_g_Yk_m == max(p_seq_g_Yk_m));
        end

    else

        % Do a partial merge based on index vector idx

        % Number of hypothesis after merge
        nm = max(idx);
        % TODO: Remove after testing
        assert(nm == length(unique(idx)));

        % Create mask to allow vectorized merging calculation
        mask = false(1, 1, nh, nm);
        ind = sub2ind(size(mask), ones(nh, 1), ones(nh, 1), (1:nh)', idx);
        mask(ind) = true;

        % Normalisation constants
        C = sum(weights .* mask, 3);

        % Merge state and output estimates
        Xk_m = reshape(sum(weights .* Xk .* mask, 3) ./ C, n, 1, []);
        if ~isempty(Yk)
            Yk_m = reshape(sum(weights .* Yk .* mask, 3) ./ C, ny, 1, []);
        else
            Yk_m = [];
        end

        % Merge state estimation error covariance
        X_diffs = Xk_m(:, :, idx) - Xk;
        to_add = pagemtimes(X_diffs, pagetranspose(X_diffs));
        Pk_m = reshape(sum(weights .* (Pk + to_add) .* mask, 3) ./ C, ...
            n, n, []);

        % TODO: May need to implement a fix here to ensure Pk_m
        % is still positive-definite. Recommended techniques
        % here: 
        %  - https://robotics.stackexchange.com/a/2004
        %
        % 1. Ensure symmetry:
        %if ~isscalar(Pk_m)
        %    Pk_m = 0.5 .* Pk_m + 0.5 .* pagetranspose(Pk_m);
        %end
        %
        % 2. Prevent underflow:
        %Pk_m = Pk_m + eps .* eye(n);
        %
        % WARNING: Before implementing these make sure everything 
        % else is working correctly!

        % Merged probabilities
        p_seq_g_Yk_m = reshape(C, [nm 1]);

        % Merged mode index
        if nargout == 5
            rk_m = nan(nm, 1);
            for i = 1:nm
                modes_to_merge = rk(idx == i);
                assert(all(modes_to_merge == modes_to_merge(1)), ...
                    "ValueError: trying to merge different modes")
                rk_m(i) = modes_to_merge(1);
            end
            % rk_m = reshape(max(reshape(double(rk), 1, 1, nh) .* mask, [], 3), nm, 1, []);
        end

    end

end