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

% Load summary results csv file
sim_label = 'popt_SP';
p_case = 1;
filename = sprintf('rod_obs_sim_%s_%d_summary.csv', sim_label, p_case);
summary_results = readtable(fullfile(results_dir, filename));

% Check what parameter values were simulated
nh_values = summary_results{:, 'MMKF_nh'};
nh_unique = unique(nh_values)';
n_min_values = summary_results{:, 'MMKF_n_min'};
n_min_unique = unique(n_min_values)';

% Check combinations are unique
combs = unique([nh_values n_min_values], 'rows');
assert(isequal(sort(combs, 1), sort([nh_values n_min_values], 1)))

fprintf("Parameter values explored:")
nh_unique
n_min_unique


%% Make plots

figure(1); clf
n_lines = numel(nh_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    nh = nh_unique(i);
    n_min_values = summary_results{nh_values == nh, 'MMKF_n_min'};
    MSE_values = summary_results{nh_values == nh, 'MSE_y_est_MMKF'};
    plot(n_min_values, MSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", nh);    
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('MSE($\hat{y}(k)$)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 650 360 180])

figure(2); clf
n_lines = numel(nh_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    nh = nh_unique(i);
    n_min_values = summary_results{nh_values == nh, 'MMKF_n_min'};
    MSE_values = summary_results{nh_values == nh, 'MSE_tr_y_est_MMKF'};
    plot(n_min_values, MSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", nh);    
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('MSE($\hat{y}(k)$) transient', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 375 360 180])


figure(3); clf
n_lines = numel(nh_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    nh = nh_unique(i);
    n_min_values = summary_results{nh_values == nh, 'MMKF_n_min'};
    MSE_values = summary_results{nh_values == nh, 'MSE_ss_y_est_MMKF'};
    plot(n_min_values, MSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", nh);    
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('MSE($\hat{y}(k)$) steady-state', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 100 360 180])