% Test pruning method

clear all;
rng(2)

% Turn this on for debugging
show_output = false;

%% Worked example

% Initialize variables
obj.n_filt = 7;
obj.n_main = 3;
obj.n_hold = 4;
obj.f_main = int16(1:obj.n_main);
obj.f_hold = int16(obj.n_main+1:obj.n_main+obj.n_hold);
obj.p_seq_g_Yk = [0.5 0.15 0.2 0.1 0 0.05 0];
nw = 2;

% Random shuffle
obj.f_main = obj.f_main(randperm(length(obj.f_main)));
obj.f_hold = obj.f_hold(randperm(length(obj.f_hold)));

if show_output
    fprintf("[%d %d %d %d][%d %d %d]\n", obj.f_hold, obj.f_main)
end

% Index of current most likely sequence
[~, f_max] = max(obj.p_seq_g_Yk);
assert(isequal(f_max, 1))

% Consistency checks - can be removed later
assert(size(obj.f_hold, 2) == obj.n_hold);
assert(size(obj.f_main, 2) == obj.n_main);
assert(isequal(sort(unique([obj.f_main obj.f_hold])), 1:obj.n_filt));

% Right-shift all filters in holding group. This causes
% the last nw values to 'roll-over' to the left of f_hold.
% e.g. circshift([1 2 3], 1) -> [3 1 2]
obj.f_hold = circshift(obj.f_hold, nw);

if show_output
    disp("After shifting holding group:")
    fprintf("[%d %d %d %d][%d %d %d]\n", obj.f_hold, obj.f_main)
end

% Filters to be moved out of holding group
f_move = obj.f_hold(1:nw);
assert(isequal(f_move, [5 4]))

% Rank hypotheses in main group according to 
% conditional probabilities
[~, i_rank] = sort(obj.p_seq_g_Yk(obj.f_main));
assert(isequal(i_rank, [1 3 2]))

% Select those with lowest probability for pruning
f_to_prune = obj.f_main(i_rank(1:nw));
assert(isequal(f_to_prune, [2 3]))

obj.f_main(i_rank(1:nw)) = f_move;

if show_output
    disp("After moving to main group:")
    fprintf("[%d %d %d %d][%d %d %d]\n", obj.f_hold, obj.f_main)
end

% Clone best sequence and put in holding group:
if show_output
    fprintf("Clone of filters{%d} -> filters{%d}, filters{%d}\n", ...
        f_max, f_to_prune)
end
obj.f_hold(1:nw) = f_to_prune;

if show_output
    disp("After replacing holding group with clone of best:")
    fprintf("[%d %d %d %d][%d %d %d]\n", obj.f_hold, obj.f_main)
end
