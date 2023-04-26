% Produce table of MKF_SF parameter optimization simulation
% results
%
% This is for the table '...' in thesis report.
%

clear all;

% Main folder where all simulation results are saved
sims_dir = 'simulations';

% Simulation group name
sim_name = "rod_obs_sim1_MKF_SF_popt";

% Create main directory for these simulations
base_dir = fullfile(sims_dir, sim_name);

% Sub-directory to store plot images
plot_dir = fullfile(base_dir, 'plots');
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Choose observer name
obs_label = "MKF_SF1";

% Specify column names for which values should be identical
id_cols = {'t_stop', 'Ts', 'nT', 'nu', 'ny', 'n', ...
    'epsilon', 'sigma_M', ...
    'sigma_wp_1_1', 'sigma_wp_1_2', ...
    'MKF_SF1_R', 'MKF_SF1_Q0_1_1', 'MKF_SF1_Q0_2_2' ...
};

% Specify column names for which values should be unique
uniq_cols = {'seed', 'MKF_SF1_nf', 'MKF_SF1_m', 'MKF_SF1_d'};

% Load simulation results summary file
filename = "rod_obs_sim_outputs_summary.csv";
results_table = readtable(fullfile(base_dir, 'results', filename));
fprintf("Existing results loaded from file: %s\n", filename)

[n_rows, n_cols] = size(results_table);
fprintf("Table size: (%d, %d)\n", n_rows, n_cols)

col_names = results_table.Properties.VariableNames;
n_obs = 0;
obs_cols = {};
n_mse = 0;
mse_cols = {};
for i = 1:numel(col_names)
    name = col_names{i};
    if startsWith(name, 'obs')
        n_obs = n_obs + 1;
        obs_cols{n_obs} = name;
    elseif startsWith(name, 'RMSE_')
        n_mse = n_mse + 1;
        mse_cols{n_mse} = name;
    end
end

% Check all simulations are with same observers and 
% other parameters match or are unique (for each seed)
% (This loop will only execute if there is more than one row)
for i = 2:n_rows

    obs_labels = results_table{i, obs_cols};
    assert(isequal(obs_labels, results_table{1, obs_cols}));

    match_params = results_table{i, id_cols};
    assert(isequal(match_params, results_table{1, id_cols}));

    uniq_params = results_table{i, uniq_cols};
    assert(~isequal(uniq_params, results_table{1, uniq_cols}));
    
end

obs_labels = unique(results_table{:, obs_cols});
fprintf("Observers: %s\n", strjoin(obs_labels));

% Iterate over each set of results
seed_values = unique(results_table.seed);
fprintf("Seeds: %s\n", strjoin(string(seed_values)));

% Select one seed
seed = 6;

% Get measurement noise value
sigma_M = results_table.sigma_M(end);


%% Find top 20 results with d > 1 (lowest y_est RMSE)

d_min = 1;

selected_cols = {
    char(obs_label + '_f'), char(obs_label + '_nf'), ...
    char(obs_label + '_m'), char(obs_label + '_d'), ...
    char(obs_label + '_nh'), char(obs_label + '_beta'), ...
    char('RMSE_y_est_' + obs_label)};

rows = (results_table.MKF_SF1_d >= d_min) & (results_table.seed == seed);
top_20 = sortrows(results_table(rows, :), ...
    'RMSE_y_est_MKF_SF1');
top_20 = top_20(1:20, selected_cols);

top_20


%% Display Latex table script

n_rows = 10;
col_labels = {'$N_f$', '$n_{max}$', '$N_d$', '$n_h$', '$\beta$', ...
    '$\text{RMSE}(\hat{Y}(N),Y(N))$'};
selected_cols = {char(obs_label + '_f'), char(obs_label + '_m'), ...
    char(obs_label + '_d'), char(obs_label + '_nh'), ...
    char(obs_label + '_beta'), char('RMSE_y_est_' + obs_label)};

fprintf("Latex table code:\n")
disp("\begin{tabular}{p{0.05\textwidth}>{\centering\arraybackslash}p{0.07\textwidth}" + ...
    ">{\centering\arraybackslash}p{0.07\textwidth}" + ...
    ">{\centering\arraybackslash}p{0.07\textwidth}" + ...
    ">{\centering\arraybackslash}p{0.07\textwidth}" + ...
    ">{\centering\arraybackslash}p{0.24\textwidth}}")
disp(strjoin(col_labels, ' & ') + "  \\")
disp("\hline")
fprintf("%% See script rod_obs_sim1_MKF_SF98_popt_table.m\n")
fprintf("%% %s results with seed = %d, sigma_M = %g, d_min = %d\n", ...
    datetime(), seed, sigma_M, d_min)
for i = 1:n_rows
    row_values = table2cell(top_20(i, selected_cols));
    fprintf("%3d & %3d & %3d & %3d & %6.4f & %6.4f \\\\\n", row_values{:})
end
disp("\hline")
disp("\end{tabular}")
		