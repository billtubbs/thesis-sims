%% Plot simulation results
%
% Author: Bill Tubbs
% 
% Input data:
%  - Simulation results in results folder
%

clear all

% Dependencies:
addpath('process-observers')
addpath('~/ml-plot-utils')

% Sub-directories used
data_dir = 'data';
results_dir = 'results';
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Specify application case
p_case = 1;  % Only case 1 used here


%% Load observers

% Load observers
% This is needed for system dimensions and model predictions
rod_obs_P1D_c4

% Select observers to include
observers = {KF1, KF2, MMKF, SKF};
n_obs = length(observers);


%% Prepare labels for tables and plots

rod_obs_make_labels


%% Make plot of process input and output data

% Select input sequence:
i_in_seq = 1;  % seq. 1 used for Fig 4 in paper.

filename = sprintf('rod_obs_sim_%d_%d.csv', p_case, i_in_seq);
opts = detectImportOptions(fullfile(results_dir, filename));
opts.VariableNamesLine = 1;
sim_data = readtable(fullfile(results_dir, filename), opts, ...
                     'ReadVariableNames', true);
fprintf("Observer simulation results loaded from file: %s\n", filename)

t = sim_data.t;
nT = size(t, 2) - 1;
if ismember('U', sim_data.Properties.VariableNames)
    U = sim_data.U;
else
    U = zeros(nT+1, 0);
end
Pd = sim_data.Pd;
Wp = sim_data.Wp;
Y = sim_data.Y;
Y_m = sim_data.Y_m;

% Calculate plant output predictions with the model
Y_model = lsim(Gpss, Wp, t);

figure(1); clf
make_iodmplot(Y, Y_m, Y_model, t, [U(:, u_meas) Pd], [u_labels p_labels], ...
    [y_labels y_m_labels y_model_labels])
xlim(t([1 end]))
set(gcf, 'Position', [100 800 360 240])
filename = sprintf('rod_obs_sim_%d_ioplot.png', i_in_seq);
saveas(gcf,fullfile(plot_dir, filename))
filename = sprintf('rod_obs_sim_%d_ioplot.pdf', i_in_seq);
saveas(gcf,fullfile(plot_dir, filename))
filename = sprintf('rod_obs_sim_%d_ioplot.svg', i_in_seq);
saveas(gcf,fullfile(plot_dir, filename))


%% Plot of observer estimates vs. true values

% Select input sequence:
i_in_seq = 3;

% Choose which observers to include in plot
obs_labels = {'KF1', 'KF2', 'MMKF'};
n_obs = numel(obs_labels);

% Load simulation data
filename = sprintf('rod_obs_sim_%d_%d.csv', p_case, i_in_seq);
sim_data = readtable(fullfile(results_dir, filename));
fprintf("Observer simulation results loaded from file: %s\n", filename)

t = sim_data.t;
nT = size(t, 1) - 1;

% Get simulation parameters from summary file
filename = sprintf('rod_obs_sim_%d_summary.csv', p_case);
summary_data = readtable(fullfile(results_dir, filename));

selection = find((summary_data.p_case == p_case) ...
    & (summary_data.i_in_seq == i_in_seq));
if numel(selection) > 1
    warning(">1 matching simulation result found. Using the last.")
end
summary_data = summary_data(selection(end), :);

% Check which observers were simulated
n_obs_sim = summary_data.n_obs;
obs_sim_labels = cell(1, n_obs_sim);
for i = 1:summary_data.n_obs
    label = sprintf("obs_%d", i);
    obs_sim_labels(i) = summary_data{1, label};
end
obs_data_exists = ismember(obs_labels, obs_sim_labels);
if ~all(obs_data_exists)
    fprintf("No data found for '%s'\n", obs_labels{~obs_data_exists})
    error("No simulation data for observer.")
end

% Get observer output estimates (ny == 1 only)
Y_est = cell(1, n_obs);
obs_sim_y_est_labels = vector_element_labels('Y_est', '', n_obs_sim);
for i = 1:n_obs
    obs_i = find(strcmp(obs_sim_labels, obs_labels{i}));
    y_est_label = obs_sim_y_est_labels{obs_i};
    Y_est{i} = sim_data{:, y_est_label};
end

% Get observer state estimates
X_est = cell(1, n_obs);
obs_sim_x_est_labels = vector_element_labels('X_est', '', n_obs_sim);
for i = 1:n_obs
    obs_i = find(strcmp(obs_sim_labels, obs_labels{i}));
    x_est_label = obs_sim_x_est_labels{obs_i};
    X_i = nan(nT+1, n);
    for j = 1:n
        x_est_label_sub = sprintf("%s_%d", x_est_label, j);
        X_i(:, j) = sim_data{:, x_est_label_sub};
    end
    X_est{i} = X_i;
end

Pd = sim_data.Pd;
%Wp = sim_data.Wp;
Y = sim_data.Y;
%Y_m = sim_data.Y_m;

figure(2); clf
axs = repmat(axes, 1, ny+nw);

for i = 1:ny
    axs(i) = subplot(ny+nw,1,i);
    labels = cell(1, n_obs);
    for j = 1:n_obs
        plot(t, Y_est{j}(:, ny),'Linewidth', 2); hold on;
        labels{j} = sprintf("$%s$ by %s", y_est_plot_labels{i}, obs_labels{j});
    end
    %stairs(t, Y_m(:, i), 'k.')
    stairs(t, Y(:, i), 'k--', 'Linewidth', 2)
    xlim(t([1 end]))
    %ylim(axes_limits_with_margin([Y Y_est], 0.1))
    ylim([-15 4])
    set(gca, 'TickLabelInterpreter', 'latex')
    %xlabel('$t$', 'Interpreter', 'Latex')
    y_label = string2latex(strjoin([y_labels(i) y_est_plot_labels(i)], ', '));
    ylabel(y_label, 'Interpreter','Latex')
    legend([obs_labels string2latex(y_labels(i))], 'Interpreter', 'Latex', 'Position', [0.39, 0.78, .13, .12])
    %legend([labels string2latex(y_labels(i))], ...
    %    'Interpreter','Latex', 'Location', 'best')
    %txt = sprintf('(%s) Observer estimates of process output %d', char(96+i), i);
    %title(txt, 'Interpreter','Latex')
    grid on
end


idx = n;  % specify states which represent input disturbances
for i = 1:nw
    axs(ny+i) = subplot(ny+nw,1,ny+i);
    labels = cell(1, n_obs);
    for j = 1:n_obs
        plot(t, X_est{j}(:, idx(i)), 'Linewidth', 2); hold on
        labels{j} = sprintf("$%s$ by %s", x_est_plot_labels{idx(i)}, obs_labels{j});
    end
    stairs(t, Pd(:, i), 'k--', 'Linewidth', 2);
    xlim(t([1 end]))
    %ylim(axes_limits_with_margin([Pd X_est(:, idx(i)+ + n*(n_obs-1))], 0.1))
    ylim([-0.2 0.5])
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel('$t$', 'Interpreter', 'Latex')
    y_label = string2latex(strjoin([p_labels(i) x_est_plot_labels(idx(i))], ', '));
    ylabel(y_label,'Interpreter','Latex')
    legend([obs_labels string2latex(p_labels(i))], 'Interpreter', 'Latex', 'Position', [0.37, 0.13, .13, .12])
    %legend([labels string2latex(p_labels(i))], 'Interpreter','Latex', 'Location', 'best')
    txt = sprintf('(%s) Observer estimates of input disturbance %d', char(96+i+2), i);
    %title(txt, 'Interpreter','Latex')
    grid on
end

linkaxes(axs, 'x')
set(gcf, 'Position', [100 600 360 360])
filename = sprintf('rod_obs_sim_%d_est.png', i_in_seq);
saveas(gcf,fullfile(plot_dir, filename))
filename = sprintf('rod_obs_sim_%d_est.pdf', i_in_seq);
saveas(gcf,fullfile(plot_dir, filename))
filename = sprintf('rod_obs_sim_%d_est.svg', i_in_seq);
saveas(gcf,fullfile(plot_dir, filename))

return


%% Plot all observer estimates over time

% TODO: Make this plot work for 2x2 systems.

% % Index to states that represent input disturbances
% idx = find(~u_meas);
% for i = 1:n_dist
% 
%     figure(2+i); clf  % Watch out if n_dist > 2
% 
%     % Plot only x2_est(k)
%     ax1 = subplot(2,1,1);
%     y_values = X_est(:, n*((1:n_obs) - 1) + idx(i));
%     plot(t, y_values); hold on
%     stairs(t, Pd(:,i), 'k-')
%     ylim(axes_limits_with_margin(y_values, 0.1, [0 1]))
%     set(gca, 'TickLabelInterpreter', 'latex')
%     xlabel('$t$', 'Interpreter', 'Latex')
%     ylabel('$p(k)$ and $\hat{x}_2(k|k-1)$', 'Interpreter', 'Latex')
%     grid on
%     labels = [obs_labels string2latex(p_labels(i))];
%     legend(labels, 'Interpreter', 'Latex', 'Location', 'best')
%     title('(a) State estimates compared to actual', 'Interpreter', 'Latex')
% 
%     ax2 = subplot(2,1,2);
%     labels = {};
%     % Calculated cumulative sum squared errors
%     y_values = cumsum((repmat(X(:, idx(i)), 1, n_obs) - y_values).^2, 1);
%     plot(t, y_values); hold on
%     max_min = [min(min(y_values)) max(max(y_values))];
%     bd = max([0.1 diff(max_min)*0.1]);
%     ylim(max_min + [-bd bd])
%     set(gca, 'TickLabelInterpreter', 'latex')
%     xlabel('$t$', 'Interpreter', 'Latex')
%     ylabel('$\sum_1^k{(x_2(k)-\hat{x}_2(k|k-1))^2}$','Interpreter','Latex')
%     grid on
%     legend(obs_labels, 'Interpreter','Latex', 'Location', 'best')
%     title('(b) Cumulative sum of squared errors', 'Interpreter', 'Latex')
% 
%     linkaxes([ax1, ax2], 'x')
% 
%     filename = sprintf('rod_obs_sim_%d_est_err.png', i_in_seq);
%     saveas(gcf,fullfile(plot_dir, filename))
% 
% end


%% Plot MKF filter results

% Plot only first input disturbance estimate
idx = 3;
i = 1;

for f = 1:n_obs_mkf

    obs = observers{observers_mkf(f)};
    % Make array of filter labels (maximum 10)
    nf_show = min(obs.n_filt, 10);
    labels = cell(1, nf_show);
    for j = 1:nf_show
        labels{j} = obs.filters{j}.label;
    end

    figure(19+f); clf

    % Plot output estimates
    ax1 = subplot(2,1,1);
    cols = ((1:nf_show) - 1)*ny + 1;
    plot(t, sim_out.MKF_Y_est{f}(:, cols)); hold on;
    stairs(t, Y(:, i), 'k-')
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel('$t$', 'Interpreter', 'Latex')
    y_label = sprintf('$%s$ and $%s$', p_labels{i}, y_est_plot_labels{1});
    ylabel(y_label, 'Interpreter', 'Latex')
    grid on
    legend([labels string2latex(y_m_labels{i})], 'Interpreter', 'Latex')
    title('(a) Output estimates', 'Interpreter', 'Latex')

    % Plot state estimates
    ax2 = subplot(2,1,2);
    cols = ((1:nf_show) - 1)*n + idx(i);
    plot(t, sim_out.MKF_X_est{f}(:, cols)); hold on;
    stairs(t, Pd(:, i), 'k-')
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel('$t$', 'Interpreter', 'Latex')
    y_label = sprintf('$%s$ and $%s$', p_labels{i}, x_est_plot_labels{idx(i)});
    ylabel(y_label, 'Interpreter', 'Latex')
    grid on
    legend([labels string2latex(p_labels{i})], 'Interpreter', 'Latex')
    title('(b) State estimates', 'Interpreter', 'Latex')

    linkaxes([ax1, ax2], 'x')

end


%% Bar plot summary of MSE estimates for each observer

figure(5)

for i = 1:ny
    subplot(ny,1,i)
    y_values = Y_mse(i:ny:end-ny+i);
    bar(categorical(obs_labels), y_values);
    set(gca, 'TickLabelInterpreter', 'latex')
    ylim(axes_limits_with_margin(y_values, 0.1, [0 0.01]))
    xlabel('Observer', 'Interpreter', 'Latex')
    ylabel('MSE', 'Interpreter', 'Latex')
    grid on
    title_txt = sprintf("(%s) Mean-squared errors of estimates of $y_%d(k)$", char(96+i), i);
    title(title_txt, 'Interpreter', 'Latex')
end

filename = sprintf('rod_obs_sim_%d_MKF_est_err_bar.png', i_in_seq);
saveas(gcf, fullfile(plot_dir, filename))

if numel(observers_mkf) == 0
    return
end


%% Plot MKF filter results - detailed plot

% Choose one MKF/AFMM filter
f = 1;

% Choose which shock
i_shock = 3;

% Identify time range to plot based on first shock occurrence
t_shocks = t(sum(alpha == 1, 2) > 0);
t_range = [t_shocks(i_shock)-Ts*20 t_shocks(i_shock)+Ts*55];
i_plot = (t >= t_range(1)) & (t <= t_range(2));

obs = observers{observers_mkf(f)};

% Make array of filter labels (only max 10 will be displayed)
nf_show = min(obs.n_filt, 10);
labels = cell(1, obs.n_filt);
for j = 1:obs.n_filt
    labels{j} = obs.filters{j}.label;
end

% Estimated states by MKF/AFMM
X_est_plot = X_est(i_plot, n*(observers_mkf(f)-1)+1:n*observers_mkf(f));
%X_err_plot = X(i_plot, :) - X_est_plot;

% Plot estimates of input disturbances after first shock

figure(7); clf
axs = repmat(axes, 1, nw);
for i = 1:nw
    axs(i) = subplot(nw,1,i);
    cols = 3:n:(obs.n_filt*n + 1);
    data = sim_out.MKF_X_est{f}(i_plot, cols(1:nf_show));
    plot(t(i_plot), data); hold on;
    plot(t(i_plot), X_est_plot(:, idx(i)), '-', 'Linewidth', 2);
    stairs(t(i_plot), Pd(i_plot,i), 'k-');
    ylim(axes_limits_with_margin(Pd(i_plot,i), 0.5))
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel('$t$', 'Interpreter', 'Latex')
    label = sprintf('$%s$ and $%s$', p_labels{i}, x_est_plot_labels{idx(i)});
    ylabel(label, 'Interpreter','Latex')
    grid on
    legend([labels(1:nf_show) {obs.label} {string2latex(p_labels{i})}], 'Interpreter', 'Latex')
    title(sprintf('(%s) Estimates of state $%s$', char(96+i), ...
        x_labels{idx(i)}), 'Interpreter', 'Latex')
end
linkaxes(axs, 'x')

filename = sprintf('rod_obs_sim_%d_%s_est_det.png', i_in_seq, obs.label);
saveas(gcf,fullfile(plot_dir, filename))

% Summary of these results
sim_results_shock = array2table( ...
    [
        t(i_plot) alpha(i_plot) ...
        round(sim_out.MKF_X_est{f}(i_plot, cols),3) ...
        round(X_est_plot(:, idx(i)), 3) ...
    ], ...
    'VariableNames', [{'t', 'alpha'} labels obs.label] ...
)


% Plot errors in state estimates

% figure(8); clf
% axs = repmat(axes, 1, n_dist);
% for i = 1:n_dist
%     axs(i) = subplot(n_dist,1,i);
%     cols = idx(i):n:(obs.n_filt*n + 1);
%     plot(t(i_plot), MKF_X_errors{f}(i_plot, cols)); hold on
%     plot(t(i_plot), X_err_plot(:, idx(i)), '-', 'Linewidth', 2);
%     set(gca, 'TickLabelInterpreter', 'latex')
%     xlabel('$t$', 'Interpreter', 'Latex')
%     label = sprintf('$%s-%s$', p_labels{i}, x_est_plot_labels{idx(i)});
%     ylabel(label,'Interpreter','Latex')
%     grid on
%     legend([labels {obs.label}], 'Interpreter', 'Latex')
%     title(sprintf('(%s) Estimatation errors of state $%s$', char(96+i), ...
%         x_labels{idx(i)}), 'Interpreter', 'Latex')
% end
% linkaxes(axs, 'x')
% 
% filename = sprintf('rod_obs_sim_%d_%s_err_det.png', i_in_seq, obs.label);
% saveas(gcf,fullfile(plot_dir, filename))


% % Plot of MSEs immediately after shocks occurred
% 
% figure(10); clf
% t_range = (0:shock_period-1)';
% axs = repmat([axes], 1, ny+n_dist);
% for i = 1:n
%     axs(i) = subplot(n,1,i);
%     for f = 1:n_obs
%         plot(t_range, X_mses_after_shocks{f}(:,i), 'Linewidth', 2);
%         hold on
%     end
%     grid on
%     ylabel('MSE','Interpreter','Latex')
%     title_txt = sprintf("(%s) Mean-squared errors of estimates of $x_%d(k)$", char(96+i), i);
%     title(title_txt, 'Interpreter', 'Latex')
%     legend(obs_labels,'Interpreter','Latex')
% end
% xlabel('$k-k_{shock}$','Interpreter','Latex')
% linkaxes(axs, 'x')
% 
% filename = sprintf('rod_obs_sim_%d_mse_ashocks.png', i_in_seq);
% saveas(gcf,fullfile(plot_dir, filename))