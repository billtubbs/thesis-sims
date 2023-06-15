% Test function indseq_from_drvs.m

assert(indseq_from_drvs(0) == 1)
assert(indseq_from_drvs(false) == 1)

assert(indseq_from_drvs(1) == 2)
assert(indseq_from_drvs(true) == 2)

assert(isequal(indseq_from_drvs([1 0 0]'), [2 1 1]))
assert(isequal(indseq_from_drvs([true false false]'), [2 1 1]))

% Combinations of 2 binary variables
G = [0     1
     0     1
     1     0
     1     1
     0     0];
assert(isequal(indseq_from_drvs(G), [3 3 2 4 1]))

% Check this matches result of shock_combinations_lte
S = cell2mat(shock_combinations_lte(2, 2));
assert(isequal( ...
    G + 1, ...
    S(indseq_from_drvs(G), :) ...
))

% All combinations of 3 binary values
G = [0 0 0
     1 0 0
     0 1 0
     0 0 1
     1 1 0
     1 0 1
     0 1 1
     1 1 1];
assert(isequal(indseq_from_drvs(G), 1:8))

% Check this matches result of shock_combinations_lte
S = cell2mat(shock_combinations_lte(3, 3));
assert(isequal( ...
    G + 1, ...
    S(indseq_from_drvs(G), :) ...
))


%% Example from docstring

% Three binary sequences
Gamma1 = [0 0 1 0 1]';
Gamma2 = [0 1 0 0 1]';
Gamma3 = [1 0 0 0 1]';

% Single shock variable
seq = indseq_from_drvs(Gamma1);
assert(isequal(seq, [1 1 2 1 2]))

% 2 shock variables
seq = indseq_from_drvs([Gamma1 Gamma2]);
assert(isequal(seq, [1 3 2 1 4]))

% Check this matches result of shock_combinations_lte
S = cell2mat(shock_combinations_lte(2, 2));
assert(isequal( ...
    [Gamma1 Gamma2] + 1, ...
    S(seq, :) ...
))

% 3 shock variables
seq = indseq_from_drvs([Gamma1 Gamma2 Gamma3]);
assert(isequal(seq, [4 3 2 1 8]))

% Check this matches result of shock_combinations_lte
S = cell2mat(shock_combinations_lte(3, 3));
assert(isequal( ...
    [Gamma1 Gamma2 Gamma3] + 1, ...
    S(seq, :) ...
))
