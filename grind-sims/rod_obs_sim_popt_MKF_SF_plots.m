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

% Dependencies:
addpath("../plot-utils")

% Sub-directories used
results_dir = 'results';
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Choose results
obs_label = "MKF_SF95";
sim_label = "popt_" + obs_label;
p_case = 1;

% Load summary results csv file
filename = sprintf("rod_obs_sim_%s_%d_summary.csv", sim_label, p_case);
summary_results = readtable(fullfile(results_dir, filename));

% Check what parameter values were simulated
f_values = summary_results{:, obs_label + "_f"};
f_unique = unique(f_values)';
m_values = summary_results{:, obs_label + "_m"};
m_unique = unique(m_values)';
d_values = summary_results{:, obs_label + "_d"};
d_unique = unique(d_values)';

% Check combinations are unique
combs = unique([f_values m_values d_values], 'rows');
assert(isequal(sort(combs, 1), sort([f_values m_values d_values], 1)))

fprintf("Parameter values explored:")
f_unique
m_unique
d_unique


%% Sort combinations from lowest to highest RMSEs

mse_label = sprintf("RMSE_y_est_%s", obs_label);

param_labels = compose(strcat(obs_label, "_%s"), ...
    ["nf", "f", "m", "d" "beta", "nm", "nh_max"]);

top_20 = sortrows(summary_results, mse_label);
top_20 = top_20(1:20, :);
disp(top_20(:, [param_labels mse_label]))


%% Make plots

figure(1); clf

n_plots = numel(m_unique);
axs = repmat(axes, 1, n_plots);

for i = 1:n_plots

    m = m_unique(i);
    m_selection = (m_values == m);
    
    axs(i) = subplot(n_plots,1,i);

    d_sel_values = unique(summary_results{m_selection, obs_label + "_d"});
    n_lines = numel(d_sel_values);
    labels = {};
    for j = 1:n_lines
        d = d_sel_values(j);
        selection = (d_values == d) & (m_values == m);
        x_data = summary_results{selection, obs_label + "_f"};
        if sum(selection) > 0
            RMSE_values = summary_results{selection, mse_label};
            plot(x_data, RMSE_values, 'o'); hold on
            labels = [labels sprintf("$d=%d$", d)];
        end
    end
    ylabel('RMSE($\hat{Y}_N,Y_N$)', 'Interpreter', 'latex');
    title(sprintf("(%s) $m=%d$", char(96+i), m), 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on
    if i == 1
        legend(labels, 'Interpreter', 'latex', 'Location', 'best');
    end
    xlabel('$f$', 'Interpreter', 'latex');

end

linkaxes(axs, 'x')
set(gcf, 'Position', [100 100 360 50+155*n_plots])

