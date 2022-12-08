%% Calculates observer metrics
%
% Author: Bill Tubbs
% 
% Calculates observer evaluation metrics summarized over
% multiple simulations.
%
% Input data:
%  - 'rod_obs_sim_1_summary.csv' in results folder
%
% To generate the simulation summary results run the
% script 'rod_obs_sim.m'
%
% NOTE: Make sure there are no duplicate results in the
% csv file. If necessary, delete the file and re-generate
% it using 'rod_obs_sim_1_summary.csv'.

clear all

% Specify path to utility functions
addpath('plot-utils/')

% Specify which simulation case
p_case = 1;  % Not currently used

% Specify which simulations to include.
% I used 6 to 15 for observer Monte Carlo simulations
i_in_seqs = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

% Choose observers to include in calculations
obs_labels = {'KF1', 'KF2', 'MKF_SF95', 'MKF_SF1', 'MKF_SP', 'SKF'};
n_obs = numel(obs_labels);

% Load simulation results from csv file
results_dir = 'results';
filename = sprintf('rod_obs_sim_%d_summary.csv', p_case);
summary_results = readtable(fullfile(results_dir, filename));
fprintf("Simulation results loaded from file: %s\n", filename)

selection = ( ...
    ismember(summary_results.i_in_seq, i_in_seqs) ...
    & summary_results.p_case == p_case ...
);
results_found = summary_results{selection, 'i_in_seq'}';
fprintf("Results for the following simulations found:\n")
disp(results_found)
if ~all(ismember(i_in_seqs, results_found))
    fprintf("No simulation results for:\n")
    disp(i_in_seqs(~ismember(i_in_seqs, results_found)))
    error("Specified simulation results missing.")
end

% Check the metrics were produced with the same assumptions
tau_ss = summary_results{selection, 'tau_ss'};
assert(all(tau_ss == tau_ss(1)))
nu = summary_results{selection, 'nu'};
assert(all(nu == nu(1)))
ny = summary_results{selection, 'ny'};
assert(all(ny == ny(1)))
ny = ny(1);

% Output variable labels
y_est_labels = vector_element_labels('y_est', '', ny, false);

% Counts of metrics
nT_Y_MSE = summary_results{selection, 'nT_Y_MSE'};
nT_Y_MSE_tr = summary_results{selection, 'nT_Y_MSE_tr'};
nT_Y_MSE_ss = summary_results{selection, 'nT_Y_MSE_ss'};

% Make summary metrics averaged over all simulations
metrics_labels = {'MSE', 'MSE_tr', 'MSE_ss', 'Var_ss', 'MSD_ss'};
n_metrics = numel(metrics_labels);
obs_metrics_labels = cell(n_metrics, n_obs * ny);
mse_table_data = nan(n_metrics, n_obs);
for f = 1:n_obs
    obs_label = obs_labels{f};

    labels = matrix_element_labels('MSE', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_MSE = summary_results{selection, labels};

    labels = matrix_element_labels('MSE_tr', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_MSE_tr = summary_results{selection, labels};

    labels = matrix_element_labels('MSE_ss', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_MSE_ss = summary_results{selection, labels};

    labels = matrix_element_labels('Var_ss', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_Var_ss = summary_results{selection, labels};
    
    labels = matrix_element_labels('MSD_ss', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_MSD_ss = summary_results{selection, labels};

    Y_MSE = sum(Y_MSE .* nT_Y_MSE) ./ sum(nT_Y_MSE);
    Y_MSE_tr = sum(Y_MSE_tr .* nT_Y_MSE_tr) ./ sum(nT_Y_MSE_tr);
    Y_MSE_ss = sum(Y_MSE_ss .* nT_Y_MSE_ss) ./ sum(nT_Y_MSE_ss);
    Y_Var_ss = sum(Y_Var_ss .* nT_Y_MSE_ss) ./ sum(nT_Y_MSE_ss);
    Y_MSD_ss = sum(Y_MSD_ss .* nT_Y_MSE_ss) ./ sum(nT_Y_MSE_ss);
    mse_table_data(:, f) = [Y_MSE; Y_MSE_tr; Y_MSE_ss; Y_Var_ss; Y_MSD_ss];
end

% Display summary table
mse_table = array2table(mse_table_data, ...
    'RowNames', {'MSE', 'MSE in transitions', 'MSE in steady-state', ...
        'Variance in steady-state', 'MSD in steady-state'}, ...
    'VariableNames', obs_labels ...
);
fprintf("Observer performance metrics\n")
disp(mse_table)

