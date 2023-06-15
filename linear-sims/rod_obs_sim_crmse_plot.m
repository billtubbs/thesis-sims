% Script to calculate convergence of RMSE over simulation duration.
clear all

addpath("../plot-utils/")

% Main folder where all simulation results are saved
sims_dir = 'simulations';

% Simulation group name
sim_name = "rod_obs_sim1_3KF_seed";

% Main simulation results directory
base_dir = fullfile(sims_dir, sim_name);

% Sub-directories with simulation outputs
results_sub_dir = 'results';
plots_sub_dir = 'plots';

% Load simulation results summary file
filename = "rod_obs_sim_outputs_summary.csv";
results_table = readtable(fullfile(base_dir, results_sub_dir, filename));
fprintf("Existing results loaded from file: %s\n", filename)

% Identify simulation results available
sim_labels = results_table.sim_label;
assert(isequal(unique(sim_labels), sim_labels))
seed_values = unique(results_table.seed);
assert(isequal(unique(seed_values), seed_values))

n_sims = numel(sim_labels);
fprintf("%d simulation results found with unique seeds:\n", n_sims)
for i = 1:n_sims
    sim_label = sim_labels{i};
    seed = results_table{strcmp(results_table.sim_label, sim_label), 'seed'};
    fprintf(" %s: %d\n", sim_label, seed);
end

% Number of samples
nT = results_table{1, 'nT'};

% Find observer labels
obs_label_cols = regexp(results_table.Properties.VariableNames, 'obs_\d*', ...
    'match', 'once');
obs_label_cols(cellfun('isempty', obs_label_cols)) = [];
obs_labels = results_table{1, obs_label_cols};
fprintf("Observers: %s\n", strjoin(obs_labels));
n_obs = numel(obs_labels);

% Prepare arrays for cumulative RMSE data
Y_crmse = cell(1, n_obs);
for i_obs = 1:n_obs
    Y_crmse{i_obs} = nan(nT+1, n_sims);
end

% Loop through each simulation and extract observer estimation errors
for i_sim = 1:n_sims
    sim_label = sim_labels{i_sim};

    % Load simulation results file  sprintf("rod_obs_sim_outputs_%s.csv")
    filename = sprintf("rod_obs_sim_outputs_%s.csv", sim_label(end-2:end));
    sim_results = readtable(fullfile(base_dir, results_sub_dir, filename));

    % Find columns with observer output errors in table
    col_names = regexp(sim_results.Properties.VariableNames, 'E_obs_\d*', ...
        'match', 'once');
    col_names(cellfun('isempty', col_names)) = [];

    % Copy desired data columns
    Y_errors = sim_results{:, col_names};

    % Calculate cumulative RMSEs for this simulation
    crmse = cum_RMSEs(Y_errors);

    % Compute convergence of MSE over simulation to check convergence
    for i_obs = 1:n_obs
        Y_crmse{i_obs}(:, i_sim) = crmse(:, i_obs);
    end

end

% Time values
t = sim_results.t;

fprintf("Final differences between min and max:\n")
for i_obs = 1:n_obs
    last_result = Y_crmse{i_obs}(end, :);
    h = max(last_result);
    l = min(last_result);
    m = median(last_result);
    fprintf("%s: RMSE = %0.4f %0.4f (%0.1f%%) / %0.4f (%0.1f%%), max - min = %.4f\n", ...
        obs_labels{i_obs}, m, l-m, (l-m)*100/m, h-m, (h-m)*100/m, h-l);
end


%% Plot time-series statistics plot

figure(1); clf
%y_labels
x_label = "Length of simulation, $N$";
y_labels = string2latex(obs_labels);
make_tsstatplot(Y_crmse, t, y_labels, x_label, nan(1, 2), 'minmax', 'median')
y_label = '$\mathrm{RMSE}(\hat{Y}(N), Y(N))$';
ylabel(y_label, 'Interpreter', 'latex')

% Standard size is [560  420]
set(gcf, 'Position', [100 730 448 336])

%filename = sprintf("%s_crmse_statsplot", sim_name);
%save_fig_to_pdf(fullfile(base_dir, plots_sub_dir, filename))
% uSed this for report
filename = sprintf("%s_crmse_statsplot.png", sim_name);
exportgraphics(gcf,fullfile(base_dir, plots_sub_dir, filename),'Resolution',300);