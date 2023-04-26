% Make plots of RMSE for different Kalman Filter tunings from
% simulation results summary file
%
% This is for the Figure 'Tuning of Kalman filter to minimise
% estimation errors' in thesis report.
%

clear all; clc

addpath('~/ml-plot-utils')

% Main folder where all simulation results are saved
sims_dir = 'simulations';

% Simulation group name
sim_name = "rod_obs_sim1_3KF_Q";

% Create main directory for these simulations
base_dir = fullfile(sims_dir, sim_name);
if ~isfolder(base_dir)
    mkdir(base_dir);
end

% Sub-directory to store plot images
plot_dir = fullfile(base_dir, 'plots');
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end


%% Load simulation summary results from csv file

% Specify column names for which values should be identical
id_cols = {'t_stop', 'Ts', 'nT', 'nu', 'ny', 'n', ...
    'epsilon', 'sigma_M', ...
    'sigma_wp_1_1', 'sigma_wp_1_2', ...
    'KF3_model_R', 'KF3_model_Q_1_1' ...
};

% Specify column names for which values should be unique
uniq_cols = {'seed', 'KF1_model_Q_2_2', 'KF2_model_Q_2_2', 'KF3_model_Q_2_2'};

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

% Find all unique seed values
seed_values = unique(results_table.seed);
fprintf("Seeds: %s\n", strjoin(string(seed_values)));
n_seeds = numel(seed_values);

% X-axis variable for all plots
x_label = 'KF3_model_Q_2_2';
x_plot_label = '$\sigma_{w_p}^2$';

% Values for x-axis
x_values = unique(results_table{:, x_label});

% Number of state variables
n = results_table{1, 'n'};

% Choose Y-axes variables for all plots
y_labels = {'RMSE_x_est_1_KF3', 'RMSE_x_est_2_KF3'};
n_plots = numel(y_labels);
y_plot_labels = {'$\mathrm{RMSE}(\hat{X}_1, X_1)$', ...
    '$\mathrm{RMSE}(\hat{X}_2, X_2)$'};
title_strings = { ...
    '(a) Root-mean-squared errors of estimates of state $x_1(k)$', ...
    '(b) Root-mean-squared errors of estimates of state $x_2(k)$'};



%% Plot RMSEs for one selected simulation

% Plot with one selected seed - this should match the seed used
% to tune the other observers.
seed = 6;

% Select RMSE data
selection = results_table.seed == seed;
selected_data = results_table(selection, [{x_label} y_labels]);
selected_data = sortrows(selected_data);
assert(isequal(selected_data{:, x_label}, x_values))

figure(10); clf

[RMSE_mins, ind] = min(selected_data{:, y_labels}, [], 1);
RMSE_mins = RMSE_mins'; ind = ind';
x_argmins = selected_data{ind, x_label};

% Display minima
fprintf("Results for seed: %d\n", seed)
disp(table(ind, x_argmins, RMSE_mins))

for i = 1:n_plots
    subplot(n_plots, 1, i);
    y_values = selected_data{:, y_labels{i}};
    semilogx(x_values, y_values, 'o-', 'Linewidth', 2, 'Markersize', 4);
    ylim([0 inf]);
    set(gca, 'TickLabelInterpreter', 'latex')
    if i == n_plots
        xlabel(x_plot_label, 'Interpreter', 'Latex')
    end
    ylabel(y_plot_labels{i}, 'Interpreter', 'Latex')
    grid on
    title(title_strings{i}, 'Interpreter', 'Latex') 
end

% Standard size is [560  420]
set(gcf, 'Position', [100 730 448 336])

filename = sprintf('%s_seed_%d', sim_name, seed);
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Plot with lines for all seeds simulated

fprintf("Best overall values (Min total RMSE(Y)):\n")
output_labels = {'RMSE_y_est_KF3'};
selected_data = results_table(:, [{x_label} output_labels]);
RMSE_sums = groupsummary(selected_data,x_label,'sum');
assert(all(RMSE_sums{:, 'GroupCount'} == n_seeds))
RMSE_sums(RMSE_sums{:, 'sum_RMSE_y_est_KF3'} == ...
    min(RMSE_sums{:, 'sum_RMSE_y_est_KF3'}), ...
    [x_label {'sum_RMSE_y_est_KF3'}])

fig = figure(11); clf

axs = cell(1, n_plots);
for i = 1:n_plots
    axs{i} = subplot(n_plots, 1, i);
end

labels = cell(1, n_seeds);
x_argmins = nan(n_seeds, numel(y_labels));
for i_seed = 1:n_seeds

    seed = seed_values(i_seed);
    selected_results = results_table(results_table.seed == seed, :);

    n = selected_results{1, 'n'};
    n_plots = numel(y_labels);
    plot_data = selected_results(:, [{x_label} y_labels]);
    plot_data = sortrows(plot_data);
    [y_mins, ind] = min(plot_data{:, y_labels});
    x_argmins(i_seed, :) = plot_data{ind, x_label}';

    x_values = plot_data{:, x_label};
    for i = 1:n_plots
        subplot(n_plots, 1, i);
        y_values = plot_data{:, y_labels{i}};
        semilogx(x_values, y_values, 'o-', 'Linewidth', 2, 'Markersize', 4);
        hold on
    end
    labels{i_seed} = sprintf("%d%", seed);
end

for i = 1:n_plots
    axes(axs{i})
    ylim([0 inf]);
    hold on
    set(gca, 'TickLabelInterpreter', 'latex')
    if j == n_plots
        xlabel(x_plot_label, 'Interpreter', 'Latex')
    end
    ylabel(y_plot_labels{i}, 'Interpreter', 'Latex')
    grid on
    title(title_strings{i}, 'Interpreter', 'Latex')
end

% leg = legend(labels, 'Interpreter', 'Latex');
% leg.Position(1) = 0.9;
% leg.Position(2) = 0.4;

% Standard size is [560  420]
set(fig, 'Position', [100 350 448 336])

filename = sprintf('%s_all_seeds', sim_name);
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Plot of averages with min-max ranges

% Select RMSE results from table
plot_data = cell(1, n_plots);
for j = 1:n_plots
    plot_data{j} = nan(numel(x_values), n_seeds);
    for i = 1:numel(x_values)
        selection = (results_table{:, x_label} == x_values(i));
        plot_data{j}(i, :) = results_table{selection, y_labels{j}};
    end
end

fig = figure(12); clf

for j = 1:n_plots

    plot_data = nan(numel(x_values), n_seeds);
    for i = 1:numel(x_values)
        selection = (results_table{:, x_label} == x_values(i));
        plot_data(i, :) = results_table{selection, y_labels{j}};
    end

    subplot(n_plots, 1, j)
    make_tsstatplot(plot_data, x_values, "", ...
        x_plot_label, nan(1, 2), 'minmax', 'median')
    set(gca,'Xscale','log')
    l = get(gca, 'Legend');
    if j < n_plots
        xlabel([])
        set(l,'visible','off')
    else
        l_pos = get(l, 'Position');
        % [0.1663 0.3538 0.3150 0.0836]
        l_pos(1) = 0.42;
        set(l, 'Position', l_pos)
    end

    ylabel(y_plot_labels{j})
    h_lines = get(gca, 'Children');
    assert(numel(h_lines) == 2)
    h_lines(1).Marker = '.';
    h_lines(1).MarkerSize = 15;

    title(title_strings{j}, 'Interpreter', 'latex')

end

% Standard size is [560  420]
set(fig, 'Position', [100 100 448 336])

filename = sprintf('%s_statplot', sim_name);
save_fig_to_pdf(fullfile(plot_dir, filename))

