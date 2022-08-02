% Test drop_empty_cols.m

clear all

test_dir = 'test_data';

if ~isfolder(test_dir)
    mkdir(test_dir)
end


%% Test 1 - normal table data

rng(0);
A = randn(5, 1);
B = randn(5, 1);
C = randi(10, 5, 2);
D = {'c1', 'c2', 'c3', 'c4', 'c5'}';
t = table(A, B, C, D);
t_drop = drop_empty_cols(t);

filename = 'test_writetable_1.csv';
writetable(t, fullfile(test_dir, filename));

filename_drop = 'test_writetable_1_drop.csv';
writetable(t_drop, fullfile(test_dir, filename_drop));

t_read = readtable(fullfile(test_dir, filename));
t_read_drop = readtable(fullfile(test_dir, filename_drop));

assert(all(abs(t{:, {'A', 'B'}} - t_read{:, {'A', 'B'}}) ...
    < 1e-14, [1 2]))
assert(isequal(t{:, 'C'}, t_read{:, {'C_1', 'C_2'}}))
assert(isequal(t{:, 'D'}, t_read{:, 'D'}))

assert(isequal(t_read, t_read_drop))


%% Test 2 - table with empty array in one column

rng(0);
A = randn(5, 1);
B = randn(5, 0);
C = randi(10, 5, 2);
D = {'c1', 'c2', 'c3', 'c4', 'c5'}';
t = table(A, B, C, D);
t_drop = drop_empty_cols(t);

assert(size(t_drop, 2) == size(t, 2) - 1)
assert(~ismember('B', t_drop.Properties.VariableNames))

filename_drop = 'test_writetable_2_fixed.csv';
writetable(t_drop, fullfile(test_dir, filename_drop));
t_read_drop = readtable(fullfile(test_dir, filename_drop));


assert(all(abs(t{:, 'A'} - t_read_drop{:, 'A'}) ...
    < 1e-14, [1 2]))
assert(~ismember('B', t_read_drop.Properties.VariableNames))
assert(isequal(t{:, 'C'}, t_read_drop{:, {'C_1', 'C_2'}}))
assert(isequal(t{:, 'D'}, t_read_drop{:, 'D'}))
