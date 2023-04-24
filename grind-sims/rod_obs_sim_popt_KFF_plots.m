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

% Dependencies:
addpath('../plot-utils')

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

RMSE_values = summary_results{:, 'RMSE_y_est_' + obs_label};
RMSE_ss_values = summary_results{:, 'RMSE_ss_y_est_' + obs_label};
RMSE_tr_values = summary_results{:, 'RMSE_tr_y_est_' + obs_label};

labels = struct();
labels.RMSE = "RMSE($\hat{Y},Y$)";
labels.RMSE_ss = "RMSE($\hat{Y},Y$) steady-state";
labels.RMSE_tr = "RMSE($\hat{Y},Y$) transient";

figure(1); clf

loglog(q_7_7_unique, RMSE_values, 'o-', 'Linewidth', 2); hold on
loglog(q_7_7_unique, RMSE_ss_values, 'o-', 'Linewidth', 2);
loglog(q_7_7_unique, RMSE_tr_values, 'o-', 'Linewidth', 2);
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('RMSE($\hat{y}(k|k)$)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf, 'Position', [100 100 360 180])
grid
legend({'Overall', 'Steady-state', 'Transitions'})


figure(2); clf

loglog(q_7_7_unique, RMSE_values, 'o-', 'Linewidth', 2);
xlabel('$n_{min}$', 'Interpreter', 'latex');
ylabel('RMSE($\hat{y}(k|k)$)', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf, 'Position', [100 350 360 180])
grid


%% Comparison with RMSEs of MKF observers

% Load RMSE results from summary results csv files
mkf_labels = ["MKF_SF95" "MKF_SF1" "MKF_SP1"];
MKF = struct();
for label = mkf_labels
    filename = sprintf("rod_obs_sim_popt_%s_1_summary.csv", label);
    results = readtable(fullfile(results_dir, filename));
    MKF.label = label;
    MKF.(label).RMSEs = results{:, 'RMSE_y_est_' + label};
    MKF.(label).RMSEs_ss = results{:, 'RMSE_ss_y_est_' + label};
    MKF.(label).RMSEs_tr = results{:, 'RMSE_tr_y_est_' + label};
end

% Find results with minimum RMSE
KF3_min = RMSE_values == min(RMSE_values);
MKF_SF95_min = MKF.MKF_SF95.RMSEs == min(MKF.MKF_SF95.RMSEs);
MKF_SF1_min = MKF.MKF_SF95.RMSEs == min(MKF.MKF_SF95.RMSEs);
MKF_SP1_min = MKF.MKF_SP1.RMSEs == min(MKF.MKF_SP1.RMSEs);

figure(3); clf

colorOrder = get(gca, 'ColorOrder');
c1 = colorOrder(1, :);
c2 = colorOrder(3, :);
c3 = colorOrder(2, :);

loglog(MKF.MKF_SF95.RMSEs_ss, MKF.MKF_SF95.RMSEs_tr, 'o', ...
    'MarkerSize', 5, 'MarkerEdgeColor', c2); hold on
loglog(MKF.MKF_SF1.RMSEs_ss, MKF.MKF_SF1.RMSEs_tr, '+', ...
    'MarkerSize', 7, 'MarkerEdgeColor', c2);
loglog(MKF.MKF_SP1.RMSEs_ss, MKF.MKF_SP1.RMSEs_tr, 'x', ...
    'MarkerSize', 7, 'MarkerEdgeColor', c2);
loglog(RMSE_ss_values, RMSE_tr_values, '.-', 'Linewidth', 2, ...
    'MarkerSize', 16, 'Color', c1);
loglog( ...
    MKF.MKF_SF95.RMSEs_ss(MKF_SF95_min), ...
    MKF.MKF_SF95.RMSEs_tr(MKF_SF95_min), ...
    'o', 'Linewidth', 2, 'MarkerSize', 7, 'MarkerEdgeColor', c3 ...
);
loglog( ...
    MKF.MKF_SF1.RMSEs_ss(MKF_SF1_min), ...
    MKF.MKF_SF1.RMSEs_tr(MKF_SF1_min), ...
    '+', 'Linewidth', 2, 'MarkerSize', 9, 'MarkerEdgeColor', c3 ...
);
loglog( ...
    MKF.MKF_SP1.RMSEs_ss(MKF_SP1_min), ...
    MKF.MKF_SP1.RMSEs_tr(MKF_SP1_min), ...
    'x', 'Linewidth', 2, 'MarkerSize', 9, 'MarkerEdgeColor', c3 ...
);
loglog( ...
    RMSE_ss_values(KF3_min), ...
    RMSE_tr_values(KF3_min), ...
    'o', 'Linewidth', 2, 'MarkerSize', 7, 'MarkerEdgeColor', c3, ...
    'MarkerFaceColor', c3 ...
);
xlim([0.9 2.5])
ylim([2 3.6])
xlabel(labels.RMSE_ss, 'Interpreter', 'latex');
ylabel(labels.RMSE_tr, 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf, 'Position', [100 600 448 336])
grid
leg_labels = [["MKF_SF95" "MKF_SF98" "MKF_SP" "KF"] ...
    compose("%s min RMSE", mkf_labels) ...
    [obs_label + " min RMSE"]];
legend(escape_latex_chars(leg_labels), 'Interpreter', 'latex')

filename = "rod_obs_sim_popt_RMSE_scatter";
save_fig_to_pdf(fullfile(plot_dir, filename))
