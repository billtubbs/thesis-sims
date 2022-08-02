% Test function axes_limits_with_margin.m

%% Vectors - without max/min ranges

lim = axes_limits_with_margin([0.71; 1.04; 0.91], 0);
assert(isequal(lim, [0.71 1.04]))

lim = axes_limits_with_margin([0.71; 1.04; 0.91], 0.2);
assert(isequal(lim, [0.71-(0.2*0.33) 1.04+(0.2*0.33)]))

lim = axes_limits_with_margin([0.71; 1.04; 0.91]);
assert(isequal(lim, [0.71-(0.1*0.33) 1.04+(0.1*0.33)]))

lim = axes_limits_with_margin([1; 1; 1]);
assert(isequal(lim, [0.9 1.1]))

lim = axes_limits_with_margin([1; 1; 1], 0);
assert(isequal(lim, [1 1]))

lim = axes_limits_with_margin([1; 1; 1], 0.1);
assert(isequal(lim, [0.9 1.1]))


%% Vectors - with min ranges

lim = axes_limits_with_margin([0.74 1.54 0.91], 0.25, [0 1]);
assert(isequal(lim, [0 1.74]))

lim = axes_limits_with_margin([0.74 1.54 0.91], 0.25, [nan 1]);
assert(isequal(lim, [0.54 1.74]))

lim = axes_limits_with_margin([0.74 1.54 0.91], 0.25, [1 2]);
assert(isequal(lim, [0.54 2]))

lim = axes_limits_with_margin([0.74 1.54 0.91], 0.25, [-2 3]);
assert(isequal(lim, [-2 3]))


%% Vectors - with max and min ranges

lim = axes_limits_with_margin([0.74; 1.54; -0.06; 0.91], 0.25, [0 1], [-2 2]);
assert(isequal(lim, [-0.46 1.94]))

lim = axes_limits_with_margin([0.74; 1.54; -0.06; 0.91], 0.25, [0 1], [-1.5 1.5]);
assert(isequal(lim, [-0.46 1.5]))

lim = axes_limits_with_margin([0.74; 1.54; -0.06; 0.91], 0.25, [nan nan], [0 1]);
assert(isequal(lim, [0 1]))



%% Test with matrix

lim = axes_limits_with_margin([-0.71 -0.71; -0.96 -1.04; -0.80 -0.91], 0.1, [0 1]);
assert(isequal(lim, [-1.0730 1]))

lim = axes_limits_with_margin([-0.71 -0.96 -0.80; -0.71 -1.04 -0.91], 0.1, [0 1]);
assert(isequal(lim, [-1.0730 1]))
