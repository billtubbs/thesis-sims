% Test group functions used by AFMM observer

clear all


%% Test function add_to_group_with_replacement

group = [1 nan nan];
new = 4;
p = [0.3 0.1 0.4 0.2];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequaln(group, [1 4 nan]))
assert(isequal(f_removed, []))

group = [nan nan nan];
new = [1 2 3];
p = [0.6 0.1 0.3];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequal(group, [1 2 3]))
assert(isequal(f_removed, []))

group = [nan nan];
new = [4 5];
p = [0.3 0.1 0.3 0.2 0.1];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequaln(group, [4 5]))
assert(isequal(f_removed, []))

group = [1 2 3];
new = 4;
p = [0.3 0.1 0.4 0.2];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequal(group, [1 4 3]))
assert(isequal(f_removed, 2))

group = [1 2 3];
new = 4;
p = [0.3 0.2 0.4 0.1];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequal(group, [1 2 3]))
assert(isequal(f_removed, 4))

group = [1 2 3];
new = 4;
p = [0.3 0.3 0.1 0.1];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert( ...
    (isequal(group, [1 2 4]) & f_removed == 3) ...
    | (isequal(group, [1 2 3]) & f_removed == 4) ...
)

group = [3 2 1];
new = [5 4];
p = [0.3 0.15 0.2 0.25 0.1];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequal(group, [3 4 1]))
assert(isequal(f_removed, [5 2]))

% Example 2 from docstring
group = [1 2 3 nan];
new = [4 5];
p = [0.4 0.1 0.25 0.05 0.2];  % 4 is lowest
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequal(group, [1 2 3 5]))
assert(isequal(f_removed, 4))

group = [nan nan 5];
new = [4 1 2];
p = [0.3 0.1 0.3 0.2 0.1];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequaln(group, [4 1 2]))
assert(isequal(f_removed, [5]))

group = [nan 2 nan];
new = [4 1 3];
p = [0.3 0.1 0.3 0.2 0.1];
[group, f_removed] = add_to_group_with_replacement(group, new, p);
assert(isequaln(group, [4 3 1]))
assert(isequal(f_removed, [2]))


%% Test function take_from_2_groups

g1 = [1 2 3];
g2 = [4 5];
[g1, g2, new] = take_from_2_groups(g1, g2, 0);
assert(isequaln(g1, [1 2 3]));
assert(isequaln(g2, [4 5]));
assert(isequal(new, []));

g1 = [1 2 3];
g2 = [4 5];
[g1, g2, new] = take_from_2_groups(g1, g2, 1);
assert(isequaln(g1, [nan 2 3]));
assert(isequaln(g2, [4 5]));
assert(isequal(new, 1));

g1 = [1 2 3];
g2 = [4 5];
[g1, g2, new] = take_from_2_groups(g1, g2, 2);
assert(isequaln(g1, [nan nan 3]));
assert(isequaln(g2, [4 5]));
assert(isequal(new, [1 2]));

g1 = [1 2 3];
g2 = [4 5];
[g1, g2, new] = take_from_2_groups(g1, g2, 3);
assert(isequaln(g1, [nan nan nan]));
assert(isequaln(g2, [4 5]));
assert(isequal(new, [1 2 3]));

g1 = [1 2 3];
g2 = [4 5];
[g1, g2, new] = take_from_2_groups(g1, g2, 4);
assert(isequaln(g1, [nan nan nan]));
assert(isequaln(g2, [nan 5]));
assert(isequal(new, [1 2 3 4]));

g1 = [1 2 3];
g2 = [4 5];
[g1, g2, new] = take_from_2_groups(g1, g2, 5);
assert(isequaln(g1, [nan nan nan]));
assert(isequaln(g2, [nan nan]));
assert(isequal(new, [1 2 3 4 5]));

g1 = [1 2 3];
g2 = [4 5];
error = false;
try
    % This should raise the following error:
    % Index exceeds the number of array elements (2).
    [g1, g2, new] = take_from_2_groups(g1, g2, 6);
catch ME
    assert(strcmp(ME.identifier, "MATLAB:badsubscript"))
    error = true;
end
assert(error)

g1 = [1 nan];
g2 = [3 nan 4 5 nan nan];
[g1, g2, new] = take_from_2_groups(g1, g2, 0);
assert(isequaln(g1, [1 nan]));
assert(isequaln(g2, [3 nan 4 5 nan nan]));
assert(isequal(new, []));

g1 = [1 nan];
g2 = [3 nan 4 5 nan nan];
[g1, g2, new] = take_from_2_groups(g1, g2, 1);
assert(isequaln(g1, [nan nan]));
assert(isequaln(g2, [3 nan 4 5 nan nan]));
assert(isequal(new, 1));

g1 = [1 nan];
g2 = [3 nan 4 5 nan nan];
[g1, g2, new] = take_from_2_groups(g1, g2, 2);
assert(isequaln(g1, [nan nan]));
assert(isequaln(g2, [nan nan 4 5 nan nan]));
assert(isequal(new, [1 3]));

g1 = [1 nan];
g2 = [3 nan 4 5 nan nan];
[g1, g2, new] = take_from_2_groups(g1, g2, 3);
assert(isequaln(g1, [nan nan]));
assert(isequaln(g2, [nan nan nan 5 nan nan]));
assert(isequal(new, [1 3 4]));

g1 = [1 nan];
g2 = [3 nan 4 5 nan nan];
[g1, g2, new] = take_from_2_groups(g1, g2, 4);
assert(isequaln(g1, [nan nan]));
assert(isequaln(g2, [nan nan nan nan nan nan]));
assert(isequal(new, [1 3 4 5]));

g1 = [1 nan];
g2 = [3 nan 4 5 nan nan];
error = false;
try
    % This should raise the following error:
    % Index exceeds the number of array elements (2).
    [g1, g2, new] = take_from_2_groups(g1, g2, 5);
catch ME
    assert(strcmp(ME.identifier, "MATLAB:badsubscript"))
    error = true;
end
assert(error)