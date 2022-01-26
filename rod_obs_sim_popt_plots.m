%% Make MMKF parameter optimization plot
%
% Make plot of MMKF MSE for different parameter values
%
% Author: Bill Tubbs
%
% Input files:
%  - rod_obs_P1D_c4.m - process model and observers
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
sim_label = 'popt';
p_case = 1;
filename = sprintf('rod_obs_sim_%s_%d_summary.csv', sim_label, p_case);
summary_results = readtable(fullfile(results_dir, filename));

% Check what parameter values were simulated
n_filt_values = summary_results{:, 'MMKF_n_filt'};
n_filt_unique = unique(n_filt_values)';
n_min_values = summary_results{:, 'MMKF_n_min'};
n_min_unique = unique(n_min_values)';

fprintf("Parameter values explored:")
n_filt_unique
n_min_unique


%% Make plots

figure(1); clf
n_lines = numel(n_filt_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    n_filt = n_filt_unique(i);
    n_min_values = summary_results{n_filt_values == n_filt, 'MMKF_n_min'};
    MSE_values = summary_results{n_filt_values == n_filt, 'MSE_y_est_MMKF'};
    plot(n_min_values, MSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", n_filt);    
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('MSE($\hat{y}(k)$)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 650 360 180])

figure(2); clf
n_lines = numel(n_filt_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    n_filt = n_filt_unique(i);
    n_min_values = summary_results{n_filt_values == n_filt, 'MMKF_n_min'};
    MSE_values = summary_results{n_filt_values == n_filt, 'MSE_tr_y_est_MMKF'};
    plot(n_min_values, MSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", n_filt);    
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('MSE($\hat{y}(k)$) transient', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 375 360 180])


figure(3); clf
n_lines = numel(n_filt_unique);
labels = cell(1, n_lines);
for i = 1:n_lines
    n_filt = n_filt_unique(i);
    n_min_values = summary_results{n_filt_values == n_filt, 'MMKF_n_min'};
    MSE_values = summary_results{n_filt_values == n_filt, 'MSE_ss_y_est_MMKF'};
    plot(n_min_values, MSE_values, '-', 'Linewidth', 2); hold on
    labels{i} = sprintf("$n_f=%d$", n_filt);    
end
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('MSE($\hat{y}(k)$) steady-state', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on
legend(labels, 'Interpreter', 'latex');
set(gcf, 'Position', [100 100 360 180])