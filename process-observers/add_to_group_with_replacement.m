function [group, f_removed] = add_to_group_with_replacement(group, new, p)
% [group, f_removed] = add_to_group_with_replacement(group, new, p)
% Adds values in vector new to vector group according to the
% following logic:
%  1. First, replace any zero values in group with values
%     from new. If all values in new were added to group,
%     then finish, returning f_replaced = [].
%  2. Rank all values in group and new combined, according
%     to values in p using: p(f) for f = [new group]. Then,
%     replace any values in group with values from new
%     if p(f_new) > p(f_group) and return f_replaced =
%     f_group.
%
%  Arguments:
%    group : vector
%        Row vector containing indices of p (no duplicates)
%        and/or zero values.
%    new : vector
%        Row vector containing only indices of p (no duplicates).
%    p : vector
%        Row vector of real values whos magnitude will be used to
%        decide which values in group to replace (if any) with
%        values in new. Length of p must be equal to total number
%        of non-zero values in group and new.
%
%  Returns:
%    group : vector
%        Updated group values (same size as before).
%    f_removed : vector
%        Values from original group and new vectors which are
%        no longer in group.
%
% Example 1:
% >> group = [1 2 3];
% >> new = 4;
% >> p = [0.3 0.1 0.4 0.2];  % 2 is lowest
% >> [group, f_removed] = add_to_group_with_replacement(group, new, p)
% 
% group =
% 
%      1     4     3
% 
% 
% f_removed =
% 
%      2
% 
% Example 2:
% >> group = [1 2 3 0];
% >> new = [4 5];
% >> p = [0.4 0.1 0.25 0.05 0.2];  % #4 is lowest
% >> [group, f_removed] = add_to_group_with_replacement(group, new, p)
% 
% group =
% 
%      1     2     3     5
% 
% 
% f_removed =
% 
%      4
% 
    n_new = numel(new);
    n_free = sum(group == 0);
    if n_free > 0
        % First fill any empty spaces (zeros) in group
        free_spaces = find(group == 0, n_new);
        group(free_spaces) = new(1:numel(free_spaces));
        new = new(n_free+1:end);
        n_new = numel(new);
    end
    if n_new > 0
        % Now, deal with any remaining values in new.
        % Rank all indices according to values in p
        combined = [group new];
        [~, f_order] = sort(p(combined));
        % Select those with lowest p for removal
        i_min = f_order(1:n_new);
        % Replaced filters returned in f_removed
        f_removed = combined(i_min);
        %[f_replace, i_replace] = ismember(group, f_removed);
        [~, i_replace] = intersect(group, f_removed, 'stable');
        group(i_replace) = setdiff(new, f_removed, 'stable');
    else
        f_removed = [];
    end
end