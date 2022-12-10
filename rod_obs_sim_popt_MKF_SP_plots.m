%% Make MKF parameter optimization plot
%
% Make plot of RMSEs of MKF observer for different parameter 
% values.
%
% Author: Bill Tubbs
%
% Input files:
%  - rod_obs_P2DcTd4.m - process model and observers
%
% Input data:
%  - Simulation outputs from Simulink model simulations.
%    See data folder
%


clear all

% Sub-directories used
results_dir = 'results';
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Choose results
obs_label = "MKF_SP1";
sim_label = "popt_" + obs_label;
p_case = 1;

% Load summary results csv file
filename = sprintf('rod_obs_sim_%s_%d_summary.csv', sim_label, p_case);
summary_results = readtable(fullfile(results_dir, filename));

% Check what parameter values were simulated
nh_values = summary_results{:, obs_label + '_nh'};
nh_unique = unique(nh_values)';
n_min_values = summary_results{:, obs_label + '_n_min'};
n_min_unique = unique(n_min_values)';

% Check combinations are unique
combs = unique([nh_values n_min_values], 'rows');
assert(isequal(sort(combs, 1), sort([nh_values n_min_values], 1)))

fprintf("Parameter values explored:")
nh_unique
n_min_unique


%% Sort combinations from lowest to highest RMSEs

rmse_label = "RMSE_y_est_" + obs_label;

param_labels = compose(strcat(obs_label, "_%s"), ...
    ["nh", "n_min"]);

top_20 = sortrows(summary_results, rmse_label);
top_20 = top_20(1:20, :);
disp(top_20(:, [param_labels rmse_label]))


%% Make plots

figure(1); clf
n_lines = numel(nh_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    nh = nh_unique(i);
    n_min_values = summary_results{nh_values == nh, obs_label + '_n_min'};
    RMSE_values = summary_results{nh_values == nh, 'RMSE_y_est_' + obs_label};
    plot(n_min_values, RMSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", nh);
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('RMSE($\hat{y}(k)$)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 650 360 180])

figure(2); clf
n_lines = numel(nh_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    nh = nh_unique(i);
    n_min_values = summary_results{nh_values == nh, obs_label + '_n_min'};
    RMSE_values = summary_results{nh_values == nh, 'RMSE_tr_y_est_' + obs_label};
    plot(n_min_values, RMSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", nh);    
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('RMSE($\hat{y}(k)$) transient', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 375 360 180])


figure(3); clf
n_lines = numel(nh_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    nh = nh_unique(i);
    n_min_values = summary_results{nh_values == nh, obs_label + '_n_min'};
    RMSE_values = summary_results{nh_values == nh, 'RMSE_ss_y_est_' + obs_label};
    plot(n_min_values, RMSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", nh);
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('RMSE($\hat{y}(k)$) steady-state', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 100 360 180])