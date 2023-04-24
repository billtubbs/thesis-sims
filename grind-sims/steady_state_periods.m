function selection = steady_state_periods(x, n, atol)
% selection = steady_state_periods(x, n, atol)
% Returns a boolean array indicating the samples in Pd that
% occur tau samples after the last change in Pd. This is useful
% for finding periods when the system is at steady-state.
%
% Arguments
%   x : vector
%       Signal data.
%   n : integer
%       Settling time in number of samples.
%   atol : absolute tolerance used to determine if a
%       change in the signal occurred.
    if nargin == 2
        atol = 0;
    end
    assert(n >= 0)
    idxs = transition_periods(x, atol);
    n_periods = numel(idxs);
    selection = true(size(x));
    for i = 1:n_periods
        idx = idxs{i};
        idx(2) = min(idx(1) + n - 1, idx(2));
        if diff(idx) >= 0
            selection(idx(1):idx(2)) = false;
        end
    end
end