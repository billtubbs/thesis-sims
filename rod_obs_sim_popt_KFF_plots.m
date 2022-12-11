%% Make KFF parameter optimization plot
%
% Make plot of RMSEs of Kalman filter for different parameter 
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
obs_label = "KF3";
sim_label = "popt_" + obs_label;
p_case = 1;

% Load summary results csv file
filename = sprintf('rod_obs_sim_%s_%d_summary.csv', sim_label, p_case);
summary_results = readtable(fullfile(results_dir, filename));

% Check what parameter values were simulated
q_7_7_values = summary_results{:, obs_label + '_model_Q_7_7'};

% Check combinations are unique
q_7_7_unique = unique(q_7_7_values);
assert(isequal(sort(q_7_7_unique, 1), sort(q_7_7_values, 1)))

fprintf("Parameter values explored:")
q_7_7_unique



%% Sort combinations from lowest to highest RMSEs

rmse_label = "RMSE_y_est_" + obs_label;

param_labels = obs_label + "_model_Q_7_7";

results_sorted = sortrows(summary_results, rmse_label);
disp(results_sorted(:, [param_labels rmse_label]))

q_min = results_sorted{1, param_labels};
fprintf("Minimum RMSE when:\nQ(7,7) = %g\n", q_min)
fprintf("sigma_wp(2) = %g\n", sqrt(q_min))


%% Make plots

figure(1); clf

RMSE_values = summary_results{:, 'RMSE_y_est_' + obs_label};
RMSE_ss_values = summary_results{:, 'RMSE_ss_y_est_' + obs_label};
RMSE_tr_values = summary_results{:, 'RMSE_tr_y_est_' + obs_label};
loglog(q_7_7_unique, RMSE_values, 'o-', 'Linewidth', 2); hold on
loglog(q_7_7_unique, RMSE_ss_values, 'o-', 'Linewidth', 2);
loglog(q_7_7_unique, RMSE_tr_values, 'o-', 'Linewidth', 2);
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('RMSE($\hat{y}(k|k)$)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf, 'Position', [100 650 360 180])
grid
legend({'Overall', 'Steady-state', 'Transitions'})


figure(2); clf

RMSE_values = summary_results{:, 'RMSE_y_est_' + obs_label};
loglog(q_7_7_unique, RMSE_values, 'o-', 'Linewidth', 2);
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('RMSE($\hat{y}(k|k)$)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf, 'Position', [100 650 360 180])
grid