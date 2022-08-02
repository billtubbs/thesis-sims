% Test objects2tablerow.m

clear all


%% Test with containers.Map

% 1 double scalars
vars = containers.Map({'a'}, {1});
T = objects2tablerow(vars);
assert(isa(T, 'table'));
assert(all(size(T) == [1 1]));
assert(T{:, 'a'} == 1);

% 2 double scalars
vars = containers.Map({'a', 'b'}, {1.1, 200});
T = objects2tablerow(vars);
assert(all(size(T) == [1 2]));
assert(T{:, 'a'} == 1.1);
assert(T{:, 'b'} == 200);

% 2 double matrices
A = ones(2,2);
B = (1:3)';
vars = containers.Map({'A', 'B'}, {A, B});
T = objects2tablerow(vars);
assert(all(size(T) == [1 numel(A) + numel(B)]));
assert(isequal(T.Properties.VariableNames, ...
    {'A_1_1', 'A_2_1', 'A_1_2', 'A_2_2', 'B_1', 'B_2', 'B_3'}));
assert(isequal(T.Variables, [reshape(A, 1, []) reshape(B, 1, [])]));

% 1 char (stored as cell array)
vars = containers.Map({'name'}, {'Fred'});
T = objects2tablerow(vars);
assert(all(size(T) == [1 1]));
assert(strcmp(T{:, 'name'}, 'Fred'));

% 1 string
vars = containers.Map({'name'}, "Fred");
T = objects2tablerow(vars);
assert(isequal(size(T), [1 1]));
assert(strcmp(T{:, 'name'}, 'Fred'));


%% Test with struct

% Simple struct with different field types
params.a = 1;
params.name = 'Fred';
params.pi = pi;
params.items = {1, 2, 3};
params.B = (1:3)';
keys = {'params'};
values = {params};
vars = containers.Map(keys, values);
T = objects2tablerow(vars);
assert(isequal(size(T), [1 9]));
assert(isequal(T{:, 'params_B_2'}, 2))
assert(isequal(T{:, 'params_items_3'}, 3))
assert(isequal(T{:, 'params_pi'}, pi))

% multi-level struct
p.x = 0.0;
p.y = 1.0;
p.z = -2.3;
params.p = p;
keys = {'params'};
values = {params};
vars = containers.Map(keys, values);
T = objects2tablerow(vars);
assert(isequal(size(T), [1 12]));
assert(isequal(T{:, 'params_p_z'}, -2.3))


%% Test with cell array

% multi-level struct
p.x = 0.0;
p.y = 1.0;
p.z = -2.3;

% cell array
keys = {'things'};
values = {{'Fred', 100, p}};
vars = containers.Map(keys, values);
T = objects2tablerow(vars);
assert(isequal(T.Properties.VariableNames, ...
    {'things_1', 'things_2', 'things_3_x', 'things_3_y', 'things_3_z'}))
assert(isequal(T{:, 'things_1'}, {'Fred'}))
assert(isequal(T{:, 'things_2'}, 100))
assert(isequal(T{:, 'things_3_z'}, -2.3))
