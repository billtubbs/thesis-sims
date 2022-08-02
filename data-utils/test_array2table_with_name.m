% Test array2table_with_name.m

X = (1:3)';
t = array2table_with_name(X.^2, 'X_sq');

assert(isequal(t.Properties.VariableNames, {'X_sq'}))
assert(isequal(t{:,:}, X.^2))

X = reshape(1:6, 3, 2);
t = array2table_with_name(X.^2, 'X_sq');

assert(isequal(t.Properties.VariableNames, {'X_sq1', 'X_sq2'}))
assert(isequal(t{:,:}, X.^2))

t = array2table_with_name(X.^2, 'X_sq', '_');
assert(isequal(t.Properties.VariableNames, {'X_sq_1', 'X_sq_2'}))
assert(isequal(t{:,:}, X.^2))