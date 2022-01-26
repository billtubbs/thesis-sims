function [indices, diffs] = transition_periods(x, atol)
% indices, diffs = transition_periods(x, atol)
% Finds the periods in a signal when it was constant
% and the differences when it changed. Assumes the 
% system was constant before the first sample period.
%
% Arguments
%   x : Column or row vector of signal data.
%   atol : Absolute tolerance used to determine if a change
%     in the signal x occurred.
%
% Returns:
%   indices : Cell array of vectors containing the start and
%     stop indices of each transition period.
%
    if nargin == 1
        atol = 0;
    end
    k_steps = find(abs(diff(x)) > atol) + 1;
    dx = diff(x);
    diffs = dx(k_steps - 1);
    n_steps = numel(k_steps);
    indices = cell(1, n_steps);
    for i = 1:n_steps
        k1 = k_steps(i);
        if i < n_steps
            k2 = k_steps(i+1) - 1;
        else
            k2 = length(x);
        end
        indices{i} = [k1 k2];
    end
end