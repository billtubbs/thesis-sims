% Test shock combinations functions
% 
%  - shock_combinations.m
%  - shock_combinations_lte.m
%  - shock_combinations_gte_prob.m
%

%% Test shock_combinations.m

S = shock_combinations(1, 0);
assert(isequal(S, {1}))

S = shock_combinations(3, 0);
assert(isequal(S, {[1 1 1]}))

S = shock_combinations(1, 1);
assert(isequal(S, {2}))

S = shock_combinations(3, 1);
S_test = [ ...
     2   1   1
     1   2   1
     1   1   2
];
assert(isequal(cell2mat(S), S_test))

% Example from docstring:
S = shock_combinations(3, 2);
S_test = [ ...
     2   2   1
     2   1   2
     1   2   2
];
assert(isequal(cell2mat(S), S_test))

S = shock_combinations(3, 3);
S_test = [ ...
     2   2   2
];
assert(isequal(cell2mat(S), S_test))


%% Test shock_combinations_lte.m

S = shock_combinations_lte(1, 0);
assert(isequal(S, {1}))

S = shock_combinations_lte(3, 0);
assert(isequal(S, {[1 1 1]}))

S = shock_combinations_lte(1, 1);
S_test = {
    1
    2
};
assert(isequal(S, S_test))

% Example from docstring:
S = shock_combinations_lte(3, 2);
S_test = [ ...
     1   1   1
     2   1   1
     1   2   1
     1   1   2
     2   2   1
     2   1   2
     1   2   2
];
assert(isequal(cell2mat(S), S_test))

S = shock_combinations_lte(3, 3);
S_test = [ ...
     1   1   1
     2   1   1
     1   2   1
     1   1   2
     2   2   1
     2   1   2
     1   2   2
     2   2   2
];
assert(isequal(cell2mat(S), S_test))


%% Test shock_combinations_gte_prob.m

% Example from docstring:
S = shock_combinations_gte_prob(3, 0.99, [0.99; 0.01]);

S_test = [ ...
   1   1   1
   2   1   1
   1   2   1
   1   1   2
];
assert(isequal(cell2mat(S), S_test))
