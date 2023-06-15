% Test rod_obs_sim_inputs.m

addpath("../process-observers/")

clear all
test_data_dir = 'test_data';


%% Test SISO system simulation setup

sys_rodin_step

obs_rodin_step

nT = 1000;
t = Ts*(0:nT)';

rng(2)
observers = {};

rod_obs_sim_inputs

test_var_names = {'alpha', 'U', 'V', 'W', 'Pd'};
assert(isequal(input_data.Properties.VariableNames, test_var_names))

% Test data filename
filename = "sim1_input_data.csv";

% Save test data if it changed
% writetable(input_data, fullfile(test_data_dir, filename))

% Compare to test data saved in CSV file
test_input_data = readtable(fullfile(test_data_dir, filename));
assert(max(abs(input_data.Variables - test_input_data.Variables), [], [1 2]) < 1e-12)


%% Test 2x2 system simulation setup

sys_rodin_step_2x2sym

obs_rodin_step_2x2

nT = 1500;
t = Ts*(0:nT)';

rng(0)
observers = {};

rod_obs_sim_inputs

test_var_names = {'alpha', 'U', 'V', 'W', 'Pd'};
assert(isequal(input_data.Properties.VariableNames, test_var_names))

% Test data filename
filename = "sim2_input_data.csv";

% Save test data if it changed
% writetable(input_data, fullfile(test_data_dir, filename))

% Compare to test data saved in CSV file
test_input_data = readtable(fullfile(test_data_dir, filename));
assert(max(abs(input_data.Variables - test_input_data.Variables), [], [1 2]) < 1e-12)
