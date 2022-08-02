function lim = axes_limits_with_margin(X, m, min_range, max_range)
% lim = axes_limits_with_margin(X, m, min_range, max_range)
% calculates new axes limits lim = [lb ub] based on the 
% extreme values in the data vector (or matrix) X but with 
% the option to add a margin m above and below that is
% proportional to (ub-lb) and to specify additional minimum 
% and maximum ranges of values which must be respected.
%
% Example:
% >> axes_limits_with_margin([0.71; 1.04; 0.91], 0.1, [0 1])
% 
% ans =
% 
%          0    1.0730
% 
    if nargin < 4
        max_range = nan(2);
    end
    if nargin < 3
        min_range = nan(2);
    end
    if nargin == 1
        m = 0.1;
    end
    [lb, ub] = bounds(X, [1 2]);
    if ub ~= lb
        margin = (ub-lb)*m;
    else
        margin = m;
    end
    lb = max(min(lb - margin, min_range(1)), max_range(1));
    ub = min(max(ub + margin, min_range(2)), max_range(2));
    lim = [lb ub];
end
