%% Plot simulation results
%
% Author: Bill Tubbs
% 
% Input data:
%  - Simulation results in results folder
%

clear all

% Dependencies:
addpath("../process-observers")
addpath("../data-utils")
addpath("../plot-utils")

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
rod_obs_P2DcTd4

% Select observers to include
observers = {KF3, MKF_SF95, MKF_SF1, MKF_SP1, SKF};
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
make_iodmplot(Y, Y_m, Y_model, t, [U(:, u_known) Pd], ...
    string2latex([u_labels p_labels]), ...
    string2latex([y_labels y_m_labels y_model_labels]), ...
    t_label)
xlim(t([1 end]))
set(gcf, 'Position', [100 800 448 336])
filename = sprintf('rod_obs_sim_%d_ioplot', i_in_seq);
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Plot of observer estimates vs. true values

% Select input sequence:
i_in_seq = 3;

% Choose which observers to include in plot
obs_labels = {'KF3', 'MKF_SF1'};
%obs_labels = {'MKF_SF95', 'MKF_SF1', 'MKF_SP1'};
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
Wp = sim_data.Wp;
Y = sim_data.Y;
Y_m = sim_data.Y_m;

figure(2); clf
axs = repmat(axes, 1, ny+nw);

for i = 1:ny
    axs(i) = subplot(ny+nw,1,i);
    labels = cell(1, n_obs);
    for j = 1:n_obs
        plot(t, Y_est{j}(:, ny),'Linewidth', 2); hold on;
        labels{j} = sprintf("$%s$ by %s", y_est_plot_labels{i}, obs_labels{j});
    end
    plot(t, Y_m(:, i), 'k.')
    plot(t, Y(:, i), 'k-')
    xlim(t([1 end]))
    %ylim(axes_limits_with_margin([Y Y_est], 0.1))
    ylim([-23 13])
    set(gca, 'TickLabelInterpreter', 'latex')
    %xlabel(t_label, 'Interpreter', 'Latex')
    y_label = string2latex(strjoin([y_m_labels(i) y_labels(i) y_est_plot_labels(i)], ', '));
    ylabel(y_label, 'Interpreter','Latex')
    legend([escape_latex_chars(obs_labels) string2latex([y_m_labels(i) y_labels(i) ])], ...
        'Interpreter', 'Latex', 'Position', [0.39, 0.78, .13, .12])
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
    stairs(t, Pd(:, i), 'k-');
    xlim(t([1 end]))
    %ylim(axes_limits_with_margin([Pd X_est(:, idx(i)+ + n*(n_obs-1))], 0.1))
    ylim([-0.25 0.6])
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel(t_label, 'Interpreter', 'Latex')
    y_label = string2latex(strjoin([p_labels(i) x_est_plot_labels(idx(i))], ', '));
    ylabel(y_label,'Interpreter','Latex')
    legend([escape_latex_chars(obs_labels) string2latex(p_labels(i))], ...
        'Interpreter', 'Latex', 'Position', [0.37, 0.13, .13, .12])
    %legend([labels string2latex(p_labels(i))], 'Interpreter','Latex', 'Location', 'best')
    txt = sprintf('(%s) Observer estimates of input disturbance %d', char(96+i+2), i);
    %title(txt, 'Interpreter','Latex')
    grid on
end

linkaxes(axs, 'x')
set(gcf, 'Position', [100 600 448 336])
filename = sprintf('rod_obs_sim_%d_est.pdf', i_in_seq);
save_fig_to_pdf(fullfile(plot_dir, filename))
