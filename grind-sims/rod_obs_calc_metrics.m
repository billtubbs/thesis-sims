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
addpath('../plot-utils')
addpath('../data-utils')

% Specify which simulation case
p_case = 1;  % Not currently used

% Specify which simulations to include.
% I used 6 to 15 for observer Monte Carlo simulations
i_in_seqs = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

% Choose observers to include in calculations
obs_labels = {'KF1', 'KF2', 'KF3', 'MKF_SF95', 'MKF_SF1', 'MKF_SP1', 'SKF'};
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
tau_ss = summary_results{selection, 'metrics_tau_ss'};
assert(all(tau_ss == tau_ss(1)))
nu = summary_results{selection, 'nu'};
assert(all(nu == nu(1)))
ny = summary_results{selection, 'ny'};
assert(all(ny == ny(1)))
ny = ny(1);

% Output variable labels
y_est_labels = vector_element_labels('y_est', '', ny, false);

% Counts of metrics
nT_Y_MSE = summary_results{selection, 'metrics_nT_Y_RMSE'};
nT_Y_MSE_tr = summary_results{selection, 'metrics_nT_Y_RMSE_tr'};
nT_Y_MSE_ss = summary_results{selection, 'metrics_nT_Y_RMSE_ss'};

% Make summary metrics averaged over all simulations
metrics_labels = {'RMSE', 'RMSE_tr', 'RMSE_ss', 'Var_ss', 'RMSD_ss'};
n_metrics = numel(metrics_labels);
obs_metrics_labels = cell(n_metrics, n_obs * ny);
rmse_table_data = nan(n_metrics, n_obs);
for f = 1:n_obs
    obs_label = obs_labels{f};

    labels = matrix_element_labels('RMSE', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_RMSE = summary_results{selection, labels};

    labels = matrix_element_labels('RMSE_tr', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_RMSE_tr = summary_results{selection, labels};

    labels = matrix_element_labels('RMSE_ss', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_RMSE_ss = summary_results{selection, labels};

    labels = matrix_element_labels('Var_ss', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_Var_ss = summary_results{selection, labels};
    
    labels = matrix_element_labels('RMSD_ss', y_est_labels, {obs_label}, '');
    labels = labels(:)';
    Y_RMSD_ss = summary_results{selection, labels};

    Y_RMSE = sum(Y_RMSE .* nT_Y_MSE) ./ sum(nT_Y_MSE);
    Y_RMSE_tr = sum(Y_RMSE_tr .* nT_Y_MSE_tr) ./ sum(nT_Y_MSE_tr);
    Y_RMSE_ss = sum(Y_RMSE_ss .* nT_Y_MSE_ss) ./ sum(nT_Y_MSE_ss);
    Y_Var_ss = sum(Y_Var_ss .* nT_Y_MSE_ss) ./ sum(nT_Y_MSE_ss);
    Y_RMSD_ss = sum(Y_RMSD_ss .* nT_Y_MSE_ss) ./ sum(nT_Y_MSE_ss);
    rmse_table_data(:, f) = [Y_RMSE; Y_RMSE_tr; Y_RMSE_ss; Y_Var_ss; Y_RMSD_ss];
end

% Display summary table
rmse_table = array2table(rmse_table_data, ...
    'RowNames', {'RMSE', 'RMSE in transitions', 'RMSE in steady-state', ...
        'Variance in steady-state', 'RMSD in steady-state'}, ...
    'VariableNames', obs_labels ...
);
fprintf("Observer performance metrics\n")
disp(rmse_table)


%% Produce Latex code

% Choose observers to include in Latex RMSE table
obs_latex = {'KF3', 'MKF_SF95', 'MKF_SF1', 'MKF_SP1', 'SKF'};

table_data = rmse_table(:, obs_latex);

% 			Metric & KF3 & MKF--SF95 & MKF--SF1 & MKF-SP & SKF \\
% 			\hline
% 			RMSE($\hat{\mathbf{Y}},\mathbf{Y}$) overall          & 11.0 & 15.9 & 3.7 & 3.5 & 2.1 \\ 
% 			RMSE($\hat{\mathbf{Y}},\mathbf{Y}$) transient       & 21.1 & 16.1 & 7.7 & 11.2 & 5.1 \\ 
% 			RMSE($\hat{\mathbf{Y}},\mathbf{Y}$) steady-state & 7.9 & 15.9 & 2.5 & 1.1 & 1.1 \\ 
% 			Var($\hat{\mathbf{Y}}$) steady-state          & 1.8 & 15.3 & 1.9 & 0.5 & 0.2 \\ 
% 			RMSD($\hat{\mathbf{Y}},\mathbf{Y}$) steady-state       & 0.0 & 16.2 & 0.5 & 0.2 & 0.0 \\ 			
% 			\hline

fprintf("Latex table code:\n")
fprintf("\\hline\n")
fprintf("%% See script rod_obs_calc_metrics.m\n")
fprintf("%% %s results with system %s, sigma_M = %g, tau_ss = %g\n", ...
    datetime(), ...
    string(summary_results{1, 'sys_name'}), ...
    summary_results{1, 'sigma_M'}, ...
    summary_results{1, 'metrics_tau_ss'})
fprintf("RMSE($\\hat{\\mathbf{Y}},\\mathbf{Y}$) overall & ")
fprintf("%s \\\\\n", strjoin(compose("%.2f", table_data{'RMSE', :}), " & "))
fprintf("RMSE($\\hat{\\mathbf{Y}},\\mathbf{Y}$) transient & ")
fprintf("%s \\\\\n", strjoin(compose("%.2f", table_data{'RMSE in transitions', :}), " & "))
fprintf("RMSE($\\hat{\\mathbf{Y}},\\mathbf{Y}$) steady-state & ")
fprintf("%s \\\\\n", strjoin(compose("%.2f", table_data{'RMSE in steady-state', :}), " & "))
fprintf("Var($\\hat{\\mathbf{Y}}$) steady-state & ")
fprintf("%s \\\\\n", strjoin(compose("%.2f", table_data{'Variance in steady-state', :}), " & "))
fprintf("RMSD($\\hat{\\mathbf{Y}},\\mathbf{Y}$) steady-state &")
fprintf("%s \\\\\n", strjoin(compose("%.2f", table_data{'RMSD in steady-state', :}), " & "))
fprintf("\\hline\n")
