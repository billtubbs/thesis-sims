function [indices, diffs] = find_step_periods(x, n1, n2, atol)
% indices = find_step_periods(x, n1, n2, atol)
% Finds any constant periods of length n2 or longer in the
% signal x, after a step change which occurred after a
% constant period of length n1 or longer. Assumes the 
% system was constant before the first sample period.
%
% Arguments:
%   x : row or column vector of signal data.
%   n1 : Minimum length of steady-state period before change
%     in number of sample periods.
%   n2 : Minimum length of period after change before the
%       next change, in number of sample periods.
%   atol : absolute tolerance used to determine if a
%       change in the signal occurred.
%
% Returns:
%   indices : cell array of beginning and end indices of step
%     periods, like {[k1 k2], [k3 k4], ... }.
%   steps : vector of step change amounts.
%
% Example
% >> idxs = find_step_periods([0 2 2 2 0 3 3 0 0 0], 2, 3);
% >> idxs{:}
% 
% ans =
% 
%      2     4
% 
% 
% ans =
% 
%      8    10
% 
% 
    if nargin == 3
        atol = 0;
    end
    assert(n1 >= 0)
    assert(n2 >= 0)
    [idxs, diffs] = transition_periods(x, atol);
    n_periods = numel(idxs);
    indices = {};
    n_prior = inf;
    keep = false(1, n_periods);
    for i = 1:n_periods
        idx = idxs{i};
        if (diff(idx) + 1 >= n2) & (n_prior >= n1)
            keep(i) = true;
        end
        n_prior = diff(idx) + 1;
    end
    diffs = diffs(keep);
    indices = idxs(keep);
end