function [g1, g2, new] = take_from_2_groups(g1, g2, n)
% [g1, g2, new] = take_from_groups(g1, g2, n)
% Takes n values from groups g1 or both g1 and g2 according
% to the following logic:
%
%  1. If there are n or more non-nan values in g1, replace
%     the first n non-non values in g1 with nans and return
%     these in new.
%  2. If there are less than n non-nan values in g1, replace
%     all the values in g1 and the first (n - numel(g1))
%     non-nan values in g2 with nans and return the combined
%     set of values in new.
%
    if n > 0
        i_g1 = find(~isnan(g1));
        i_remove = i_g1(1:min(n, numel(i_g1)));
        new = g1(i_remove);
        g1(i_remove) = nan;
        if numel(new) < n
            i_g2 = find(~isnan(g2));
            i_remove = i_g2(1:n-numel(new));
            new = [new g2(i_remove)];
            g2(i_remove) = nan;
        end
    else
        new = [];
    end
