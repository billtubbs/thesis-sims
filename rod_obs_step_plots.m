%% Make step response plots
%
% Make min-max and mean plots of stochastic responses of
% observers to step disturbances.
%
% Author: Bill Tubbs
%
% Input data:
%  - Simulation outputs from Simulink model simulations.
%    See data folder
%

clear all
rng(0)

% Dependencies:
addpath('process-observers')
addpath('data-utils')
addpath('plot-utils')

% Sub-directories used
data_dir = 'data';
results_dir = 'results';
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Specify application case
p_case = 1;  % Only case 1 used here

% Choose sequences from which to extract step responses
i_in_seqs = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

% Specify length of responses and sample period
nT_resp = 51;
Ts = 3/60;  % hours
t_resp = Ts*(0:nT_resp-1)';

% Choose how many responses to include in plot
n_resp = 100;

% Make labels
ny = 1;
y_labels = vector_element_labels('y', '(k)', ny, false);

% Get simulation parameters from summary file
filename = sprintf('rod_obs_sim_%d_summary.csv', p_case);
summary_data = readtable(fullfile(results_dir, filename));

selection = find((summary_data.p_case == p_case) & ...
    ismember(summary_data.i_in_seq, i_in_seqs));
n_seq = numel(selection);
fprintf("%d simulation results match criteria.\n", n_seq)

% Drop any results that don't match
summary_data = summary_data(selection, :);

% Check which observers were simulated
assert(all(summary_data.n_obs == summary_data.n_obs(1)))
n_obs_sim = summary_data.n_obs(1);
obs_sim_labels = cell(1, n_obs_sim);
for i = 1:summary_data.n_obs
    label = sprintf("obs_%d", i);
    obs_sim_labels(i) = summary_data{1, label};
end

% Select observers to include
obs_labels = {'KF3', 'MKF_SF95'};
n_obs = numel(obs_labels);

obs_data_exists = ismember(obs_labels, obs_sim_labels);
if ~all(obs_data_exists)
    fprintf("No data found for '%s'\n", obs_labels{~obs_data_exists})
    error("No simulation data for observer.")
end
obs_sel = cell(1, n_obs);
for f = 1:n_obs
    if sum(strcmp(obs_labels{f}, obs_labels)) == 1
        obs_sel{f} = true;
    else
        obs_sel{f} = false;
    end
end

% Load existing results and combine
filename = sprintf('rod_obs_sim_%d_resps.csv', p_case);
resp_data = readtable(fullfile(results_dir, filename));
n_resp_sim = size(resp_data, 1);
fprintf("%d step responses loaded from file: %s\n", n_resp_sim, filename)

% Cut down to required number
selection = datasample(1:n_resp_sim, n_resp, 'Replace', false);
resp_data = resp_data(selection, :);
fprintf("%d step responses selected.\n", size(resp_data, 1))

% Size of steps
diffs = resp_data{:, 'diffs'};

% First find positive step trajectories
selection = diffs >= 0;
y_resp_labels = vector_element_labels('yk', '', nT_resp);
Y_resp = resp_data{selection, y_resp_labels}';
Y_med = nanmedian(Y_resp, 2);

% Responses of selected observers
Y_est_resp = cell(1, n_obs);
Y_est_med = cell(1, n_obs);
y_est_labels_obs = cell(1, n_obs);
f2 = 1;
for f = 1:n_obs
    if obs_sel{f} == false
        continue
    end
    label = strcat('yk_est_', obs_labels{f});
    Y_est_resp_labels = vector_element_labels(label, '', nT_resp);
    Y_est_resp{f2} = resp_data{selection, Y_est_resp_labels}';
    Y_est_med{f2} = nanmedian(Y_est_resp{f2}, 2);
    y_est_labels_obs{f2} = ['$\hat{y}(k)$ ' escape_latex_chars(obs_labels{f})];
    f2 = f2 + 1;
end

figure(3); clf
f = 1;
plot(t_resp, Y_resp, 'k-'); hold on
plot(t_resp, Y_med, 'k-', 'Linewidth', 1);
plot(t_resp, Y_est_resp{f}, 'r-');
plot(t_resp, Y_est_med{f}, 'r-', 'Linewidth', 2);
%plot(t_resp, reshape(Y_est_resps(:, 2, :), 29, 12), 'r-');
%plot(t_resp, reshape(Y_est_resps(:, 3, :), 29, 12), 'b-');
grid on

figure(4); clf
labels = [string2latex(y_labels) y_est_labels_obs];
make_tsstatplot([{Y_resp} Y_est_resp], t_resp, labels, '$t-t_{step}$');
set(gcf, 'Position', [100 400 448 336])
filename = "rod_obs_sim_resp_plot1.pdf";
save_fig_to_pdf(fullfile(plot_dir, filename))

% Plot all negative step trajectories
selection = diffs < 0;
y_resp_labels = vector_element_labels('yk', '', nT_resp);
Y_resp = resp_data{selection, y_resp_labels}';
Y_med = nanmedian(Y_resp, 2);

% Responses of selected observers
Y_est_resp = cell(1, n_obs);
Y_est_med = cell(1, n_obs);
y_est_labels_obs = cell(1, n_obs);
f2 = 1;
for f = 1:n_obs
    if obs_sel{f} == false
        continue
    end
    label = strcat('yk_est_', obs_labels{f});
    Y_est_resp_labels = vector_element_labels(label, '', nT_resp);
    Y_est_resp{f2} = resp_data{selection, Y_est_resp_labels}';
    Y_est_med{f2} = nanmedian(Y_est_resp{f2}, 2);
    y_est_labels_obs{f2} = ['$\hat{y}(k)$ ' escape_latex_chars(obs_labels{f})];
    f2 = f2 + 1;
end

figure(5); clf
f = 1;
t_resp = Ts*(0:nT_resp-1)';
plot(t_resp, Y_resp, 'k-'); hold on
plot(t_resp, Y_med, 'k-', 'Linewidth', 1);
plot(t_resp, Y_est_resp{f}, 'r-');
plot(t_resp, Y_est_med{f}, 'r-', 'Linewidth', 2);
%plot(t_resp, reshape(Y_est_resps(:, 2, :), 29, 12), 'r-');
%plot(t_resp, reshape(Y_est_resps(:, 3, :), 29, 12), 'b-');
grid on

figure(6); clf
labels = [string2latex(y_labels) y_est_labels_obs];
make_tsstatplot([{Y_resp} Y_est_resp], t_resp, labels, '$t-t_{step}$');
set(gcf, 'Position', [100 200 448 336])
filename = "rod_obs_sim_resp_plot2.pdf";
save_fig_to_pdf(fullfile(plot_dir, filename))




