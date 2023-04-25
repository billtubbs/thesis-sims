%% Plot simulation results
% 
% Dependencies:
%   - Run this after generating the simulation results in the
%     specified directory.  E.g.:
% 
%   - results/rod_obs_sim1_all_seed
%

clear all

addpath('~/ml-data-utils')
addpath('~/ml-plot-utils')

% Main folder where all simulation results are saved
sims_dir = 'simulations';

% Simulation group names
%sim_name = "rod_obs_sim1_all_seed";
sim_name = "rod_obs_sim2_all_seed";

% Choose seed to select which simulation results to plot
switch sim_name
    case "rod_obs_sim1_all_seed"
        seed = 6;  % 0 or 6 for SISO system
    case "rod_obs_sim2_all_seed"
        seed = 6;  % 0 for 2x2 system #2
end


%% Load summary simulation results

% Sub-directory where simulation results are stored
results_sub_dir = 'results';

% Main directory for main simulations
base_dir = fullfile(sims_dir, sim_name);
filename = 'rod_obs_sim_outputs_summary.csv';
fprintf("Loading summary results from '%s'...\n", base_dir)
filepath = fullfile(base_dir, results_sub_dir, filename);
% Getting a warning but can't figure out how to avoid it.
%opts = detectImportOptions(filepath);
%opts = setvaropts(opts,'Time','InputFormat','uuuu-MM-dd HH:mm:ss'); 
%results_summary = readtable(filepath, opts);
results_summary = readtable(filepath);

% Main directory for 'no noise' simulations
results_sub_dir_no_noise = 'results_no_noise';
filename = 'rod_obs_sim_outputs_summary.csv';
filepath = fullfile(base_dir, ...
    results_sub_dir_no_noise, filename);
fprintf("Loading 'no-noise' summary results from '%s'...\n", base_dir)
results_summary_no_noise = readtable(filepath);

% Check main and 'no noise' simulations match
assert(isequal(sort(results_summary.Properties.VariableNames), ...
    sort(results_summary_no_noise.Properties.VariableNames)))
selected_row = results_summary.seed == seed;
assert(sum(selected_row) == 1)
assert(results_summary{selected_row, "t_stop"} == results_summary_no_noise.t_stop)
assert(results_summary{selected_row, "Ts"} == results_summary_no_noise.Ts)
assert(results_summary{selected_row, "n"} == results_summary_no_noise.n)
sys_params_cols = startsWith(results_summary.Properties.VariableNames, "params");
assert(isequal(results_summary{selected_row, sys_params_cols}, results_summary_no_noise{:, sys_params_cols}))

% Sub-directory to store plot images
plot_dir = fullfile(base_dir, 'plots');
if ~isfolder(plot_dir)
    mkdir(plot_dir);
end

% Run checks for uniqueness, equal vars etc.
assert(isequal(unique(results_summary.seed), sort(results_summary.seed)))

col_names = results_summary.Properties.VariableNames;
n_obs = 0;
obs_cols = {};
Y_est_cols = {};
for i = 1:numel(col_names)
    name = col_names{i};
    if startsWith(name, 'obs')
        n_obs = n_obs + 1;
        obs_cols{n_obs} = name;
    end
end
obs_labels = results_summary{1, obs_cols};

% Check observer parameters too.
for label = obs_labels
    if startsWith(label, "KF")
        n = results_summary{selected_row, "n"};
        switch n
            case 1
                col = sprintf("%s_model_Q", label{1});
                assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
                col = sprintf("%s_P0", label{1});
                assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
            otherwise
                for i = 1:n
                    col = sprintf("%s_model_Q_%d_%d", label{1}, i, i);
                    assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
                    col = sprintf("%s_P0_%d_%d", label{1}, i, i);
                    assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
                end
        end
        ny = results_summary{selected_row, "ny"};
        switch ny
            case 1
                col = sprintf("%s_model_R", label{1});
                assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
            otherwise
                for i = 1:ny
                    col = sprintf("%s_model_R_%d_%d", label{1}, i, i);
                    assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
                end
        end
    elseif startsWith(label, "MKF")
        col = label + "_nh";
        assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
        n = results_summary{selected_row, "n"};
        for i = 1:n
            col = sprintf("%s_Q0_%d_%d", label{1}, n, n);
            assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
            col = sprintf("%s_P0_%d_%d", label{1}, n, n);
            assert(results_summary{selected_row, col} == results_summary_no_noise{1, col})
        end
    end
end


%% Load main set of simulation results

i_sim = find(results_summary.seed == seed);
selected_sim = results_summary(i_sim, :);
assert(size(selected_sim, 1) == 1)
sim_label = selected_sim.sim_label;
nT = selected_sim{1, 'nT'};
Ts = selected_sim{1, 'Ts'};
nu = selected_sim.nu;
n = selected_sim.n;
ny = selected_sim.ny;
u_known = boolean(selected_sim{:, ...
    col_names(startsWith(col_names, 'u_known'))});
nw = sum(~u_known);
nu = nu - nw;

% This is a bit of a hack
items = split(sim_label,'_');
sim_num = str2double(items(end));
filename = sprintf("rod_obs_sim_outputs_%03d.csv", sim_num);

fprintf("Loading simulation results from '%s'...\n", filename)
sim_results = readtable(fullfile(base_dir, results_sub_dir, filename));
sim_col_names = sim_results.Properties.VariableNames;

t = sim_results{:, 't'};
fprintf("  Selected seed: %g\n", seed)


%% Prepare plot labels

rod_obs_sim_labels


%% Plot of process inputs and outputs

% Select time range
t_stop = 300;
t_range = sim_results{t <= t_stop, 't'};

switch ny
    case 1
        Y = sim_results{t <= t_stop, 'Y'};
        Y_m = sim_results{t <= t_stop, 'Y_m'};
    otherwise
        Y = sim_results{t <= t_stop, compose('Y_%d', 1:ny)};
        Y_m = sim_results{t <= t_stop, compose('Y_m_%d', 1:ny)};
end

U = sim_results{t <= t_stop, sim_col_names(startsWith(sim_col_names, 'U'))};
Pd = sim_results{t <= t_stop, sim_col_names(startsWith(sim_col_names, 'Pd'))};
alpha = sim_results{t <= t_stop, sim_col_names(startsWith(sim_col_names, 'alpha'))};
Wp = U(:, ~u_known);

figure(1); clf
make_iodplot(Y, Y_m, t_range, [U(:, u_known) Pd], ...
    string2latex([u_labels p_labels]), string2latex([y_labels yM_labels]), ...
    time_label)

% Standard size is [560  420]
set(gcf, 'Position', [50 950 448 336])

filename = sprintf('%s_ioplot', sim_name);
save_fig_to_pdf(fullfile(plot_dir, filename))

%% Choose observers to include in plots

% Select observers
%obs_sel_labels = {'KF1', 'KF2', 'KF3'};   % for sim1
%obs_sel_labels = {'KF3', 'SKF', 'MKF_SF95'};
%obs_sel_labels = {'KF3', 'SKF', 'MKF_SP1'};
%obs_sel_labels = {'KF1', 'KF3'};
obs_sel_labels = {'MKF_SF95', 'MKF_SF1', 'SKF'};
%obs_sel_labels = {'MKF_SP1', 'SKF'};
%obs_sel_labels = {'KF3', 'MKF_SF95', 'MKF_SP1', 'SKF'};
%obs_sel_labels = {'KF3', 'MKF_SF95', 'MKF_SF1', 'MKF_SP1', 'SKF'};  % for sim1 & sim3

obs_sel = find(ismember(obs_labels, obs_sel_labels));
n_obs = numel(obs_sel);


%% Plot of observer estimates vs. true values

% Get observer estimates
X_est = sim_results{t <= t_stop, ...
    startsWith(sim_results.Properties.VariableNames, 'X_est')};
Y_est = sim_results{t <= t_stop, ...
    startsWith(sim_results.Properties.VariableNames, 'Y_est')};

% Measurements and true signals
Y_m = sim_results{t <= t_stop, sim_col_names(startsWith(sim_col_names, 'Y_m'))};
Y = sim_results{t <= t_stop, sim_col_names(startsWith(sim_col_names, 'Y'))};
Pd = sim_results{t <= t_stop, sim_col_names(startsWith(sim_col_names, 'Pd'))};

figure(3); clf
axs = repmat(axes, 1, ny+nw);

for i = 1:ny
    axs(i) = subplot(ny+nw,1,i);
    cols = i + (obs_sel-1)*ny;
    y_values = Y_est(:, cols);
    labels = {};
    for j = obs_sel
        plot(t_range, Y_est(:, ny*(j-1)+i), 'Linewidth', 2); hold on;
        label = sprintf("$%s$ by %s", y_est_plot_labels{i}, ...
            escape_latex_chars(obs_labels{j}));
        labels = [labels {label}];
    end
    %plot(t_range, Y_m(:, i), 'k.')
    plot(t_range, Y(:, i), 'k-')
    ylim(axes_limits_with_margin([Y_m(:,i) Y(:,i) y_values], 0.1))
    set(gca, 'TickLabelInterpreter', 'latex')
    %xlabel(time_label, 'Interpreter', 'Latex')
    y_label = string2latex(strjoin([y_est_plot_labels(i) y_labels(i)], ', '));
    ylabel(y_label, 'Interpreter','Latex')
    if i == 1  % legend on first plot only
        leg = legend([labels string2latex(y_labels(i))], ...
            'Interpreter', 'Latex', 'Location', 'best');
    end
    %set(leg, 'Position', [0.415 0.762 0.2624 0.1449])
    txt = sprintf('(%s) Observer estimates of process output', char(96+i));
    if ny > 1
        txt = strjoin({txt, char(string(i))});
    end
    title(txt, 'Interpreter','Latex')
    grid on
end

idx = find(~u_known);
for i = 1:nw
    axs(ny+i) = subplot(ny+nw,1,ny+i);
    cols = idx(i) + (obs_sel-1)*n;
    y_values = X_est(:, cols);
    labels = {};
    y_lims = axes_limits_with_margin(y_values(10:end, :), 0.1);  % avoid large initial errors
    for j = obs_sel
        plot(t_range, X_est(:, idx(i) + n*(j-1)), 'Linewidth', 2); hold on
        label = sprintf("$%s$ by %s", x_est_plot_labels{idx(i)}, ...
            escape_latex_chars(obs_labels{j}));
        labels = [labels {label}];
    end
    stairs(t_range, Pd(:, i), 'k-'); hold on;
    ylim(axes_limits_with_margin(y_values(10:end, :), 0.1))  % avoid large initial errors
    set(gca, 'TickLabelInterpreter', 'latex')
    if i == nw
        xlabel(time_label, 'Interpreter', 'Latex')
    end
    y_label = string2latex(strjoin([p_labels(i) x_est_plot_labels(idx(i))], ', '));
    ylabel(y_label,'Interpreter','Latex')
    if i == 1  % legend on first plot only
        leg = legend([labels string2latex(p_labels(i))], ...
            'Interpreter', 'Latex', 'Location', 'best');
    end
    %set(leg, 'Position', [0.415 0.29 0.2624 0.1449])
    txt = sprintf('(%s) Observer estimates of input disturbance %d', char(96+i+ny));
    if nw > 1
        txt = strjoin({txt, char(string(i))});
    end
    title(txt, 'Interpreter','Latex')
    grid on
end

linkaxes(axs, 'x')

% Standard size is [560  420]
set(gcf, 'Position', [50 100 448 168*(ny+nw)])

filename = sprintf('%s_y_est1', sim_name);
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Load no-noise simulation results

% no_noise = struct();
% i_sim = find(results_summary_no_noise.seed == seed);
% no_noise.selected_sim = results_summary_no_noise(i_sim, :);
% assert(size(selected_sim, 1) == 1)
% no_noise.sim_label = no_noise.selected_sim.sim_label;
% no_noise.nT = selected_sim{1, 'nT'};
% no_noise.Ts = selected_sim{1, 'Ts'};
% no_noise.nu = selected_sim.nu;
% no_noise.n = selected_sim.n;
% no_noise.ny = selected_sim.ny;
% no_noise.u_known = boolean(no_noise.selected_sim{:, ...
%     col_names(startsWith(col_names, 'u_known'))});
% no_noise.nw = sum(~u_known);
% no_noise.nu = nu - nw;
% 
% % This is a bit of a hack
% items = split(sim_label,'_');
% sim_num = str2double(items(end));
% filename = sprintf("rod_obs_sim_outputs_%03d.csv", sim_num);
% 
% fprintf("Loading 'no noise' simulation results from '%s'...\n", ...
%     filename)
% no_noise.sim_results = readtable(fullfile(base_dir, results_sub_dir_no_noise, filename));
% no_noise.sim_col_names = sim_results.Properties.VariableNames;
% 
% no_noise.t = no_noise.sim_results{:, 't'};
% 
% % Check simulation inputs are the same(those that should be)
% switch sim_name
%     case "rod_obs_sim1_all_seed"
%         assert(isequal(sim_results.alpha, no_noise.sim_results.alpha))
%         assert(isequal(sim_results.Pd, no_noise.sim_results.Pd))
%         assert(all(no_noise.sim_results{:, 'V'} == 0))
%     case {"rod_obs_sim2_all_seed", "rod_obs_sim3_all_seed"}
%         assert(isequal(sim_results.alpha_1, no_noise.sim_results.alpha_1))
%         assert(isequal(sim_results.alpha_2, no_noise.sim_results.alpha_2))
%         assert(isequal(sim_results.Pd_1, no_noise.sim_results.Pd_1))
%         assert(isequal(sim_results.Pd_2, no_noise.sim_results.Pd_2))
%         assert(all(no_noise.sim_results{:, 'V_1'} == 0))
%         assert(all(no_noise.sim_results{:, 'V_1'} == 0))
%     otherwise
%         error("Unrecognized sim_name")
% end


%% Plot of observer estimates with and without noise
 
% no_noise.X_est = no_noise.sim_results{t <= t_stop, ...
%     startsWith(no_noise.sim_results.Properties.VariableNames, 'X_est')};
% no_noise.Y_est = no_noise.sim_results{t <= t_stop, ...
%     startsWith(no_noise.sim_results.Properties.VariableNames, 'Y_est')};
% no_noise.Y_m = no_noise.sim_results{t <= t_stop, sim_col_names(startsWith(sim_col_names, 'Y_m'))};
% no_noise.Pd = no_noise.sim_results{t <= t_stop, sim_col_names(startsWith(sim_col_names, 'Pd'))};
% 
% for j = obs_sel
% 
%     figure(9+j); clf
%     axs = repmat(axes, 1, ny+nw);
% 
%     for i = 1:ny
%         axs(i) = subplot(ny+nw,1,i);
%         labels = {};
%         plot(t_range, Y_est(:, ny*(j-1)+i), 'Linewidth', 2); hold on
%         plot(t_range, no_noise.Y_est(:, ny*(j-1)+i), 'Linewidth', 2)
%         label1 = sprintf("$%s$ with noise", y_est_plot_labels{i});
%         label2 = sprintf("$%s$ no noise", y_est_plot_labels{i});
%         labels = [labels {label1, label2}];
%     
%         plot(t_range, Y_m(:, i), 'k.')
%         ylim(axes_limits_with_margin([Y_m(:,i) Y_est(:,i)], 0.1))
%         set(gca, 'TickLabelInterpreter', 'latex')
%         %xlabel(time_label, 'Interpreter', 'Latex')
%         y_label = string2latex(strjoin([yM_labels(i) y_est_plot_labels(i)], ', '));
%         ylabel(y_label, 'Interpreter','Latex')
%         if i == 1  % legend on first plot only
%             leg = legend([labels string2latex(yM_labels(i))], ...
%                 'Interpreter', 'Latex', 'Location', 'best');
%         end
%         %set(leg, 'Position', [0.415 0.762 0.2624 0.1449])
%         txt = sprintf('(%s) %s estimates of process output', char(96+i), ...
%             escape_latex_chars(obs_labels{j}));
%         if ny > 1
%             txt = strjoin({txt, char(string(i))});
%         end
%         title(txt, 'Interpreter','Latex')
%         grid on
%     end
% 
%     idx = find(~u_known);
%     for i = 1:nw
%         axs(ny+i) = subplot(ny+nw,1,ny+i);
%         labels = {};
%         y_values1 = X_est(:, idx(i) + n*(j-1));
%         y_values2 = no_noise.X_est(:, idx(i) + n*(j-1));
%         plot(t_range, y_values1, 'Linewidth', 2); hold on
%         plot(t_range, y_values2, 'Linewidth', 2)
%         label1 = sprintf("$%s$ with noise", x_est_plot_labels{idx(i)});
%         label2 = sprintf("$%s$ no noise", x_est_plot_labels{idx(i)});
%         labels = [labels {label1, label2}];
%         stairs(t_range, Pd(:, i), 'k-'); hold on;
%         ylim(axes_limits_with_margin([y_values1(10:end, :) y_values2(10:end, :)], 0.1))  % avoid large initial errors
%         set(gca, 'TickLabelInterpreter', 'latex')
%         if i == nw
%             xlabel(time_label, 'Interpreter', 'Latex')
%         end
%         y_label = string2latex(strjoin([p_labels(i) x_est_plot_labels(idx(i))], ', '));
%         ylabel(y_label,'Interpreter','Latex')
%         if i == 1  % legend on first plot only
%             leg = legend([labels string2latex(p_labels(i))], ...
%                 'Interpreter', 'Latex', 'Location', 'best');
%         end
%         %set(leg, 'Position', [0.415 0.29 0.2624 0.1449])
%         txt = sprintf('(%s) %s estimates of input disturbance %d', ...
%             char(96+i+ny), escape_latex_chars(obs_labels{j}));
%         if nw > 1
%             txt = strjoin({txt, char(string(i))});
%         end
%         title(txt, 'Interpreter','Latex')
%         grid on
%     end
%     
%     linkaxes(axs, 'x')
%     
%     % Standard size is [560  420]
%     set(gcf, 'Position', [50 615-j*75 448 168*(ny+nw)])
%     
%     filename = sprintf('%s_y_est_wout_noise_%s', sim_name, obs_labels{j});
%     save_fig_to_pdf(fullfile(plot_dir, filename))
% 
% end


%% Custom plot of observer estimates vs. true values

% Plot observer estimates of disturbances during a
% selected period

% Select specific plot ranges for each sim case
switch sim_name
    case {"rod_obs_sim1_all_seed", "rod_obs_sim1_no_noise"}
        t_sel = (t >= 90) & (t <= 120);
        y_lims = {[-0.6 1.5], [-1.6 0.5]};
        leg_loc = 'best';
    case {"rod_obs_sim2_all_seed", "rod_obs_sim2_no_noise"}
        t_sel = (t >= 115) & (t <= 175);
        y_lims = {[0.5 5.2], [-5.2 -1.6], [-1 3.7], [-3.8 -1.8]};
        leg_loc = 'best';
    case {"rod_obs_sim3_all_seed", "rod_obs_sim3_no_noise"}
        t_sel = (t >= 260) & (t <= 350);
        y_lims = {[-2.5 0.8], [-1.2 1.1], [-2 1.2], [-1.2 1.2]};
        leg_loc = 'best';
end

X_est = sim_results{t_sel, ...
    startsWith(sim_results.Properties.VariableNames, 'X_est')};
Y_est = sim_results{t_sel, ...
    startsWith(sim_results.Properties.VariableNames, 'Y_est')};
Y = sim_results{t_sel, sim_col_names(startsWith(sim_col_names, 'Y'))};
Y_m = sim_results{t_sel, sim_col_names(startsWith(sim_col_names, 'Y_m'))};
U = sim_results{t_sel, sim_col_names(startsWith(sim_col_names, 'U'))};
Pd = sim_results{t_sel, sim_col_names(startsWith(sim_col_names, 'Pd'))};
t_range2 = t(t_sel);

figure(4); clf
axs = repmat(axes, 1, ny+nw);

i_plot = 1;
for i = 1:ny
    axs(i) = subplot(ny+nw,1,i);
    labels = {};
    for j = obs_sel
        plot(t_range2, Y_est(:, ny*(j-1)+i), 'Linewidth', 2); hold on;
        label = sprintf("$%s$ by %s", y_est_plot_labels{i}, ...
            escape_latex_chars(obs_labels{j}));
        labels = [labels {label}];
    end
    stairs(t_range2, Y_m(:, i), 'k.')
    plot(t_range2, Y(:, i), 'k-')  % optional
    xlim(t_range2([1 end])')
    ylim(y_lims{i_plot}); i_plot = i_plot + 1;
    set(gca, 'TickLabelInterpreter', 'latex')
    %xlabel(time_label, 'Interpreter', 'Latex')
    y_label = string2latex(strjoin([y_est_plot_labels(i) yM_labels(i) y_labels(i)], ', '));
    ylabel(y_label, 'Interpreter','Latex')
    if i == 1  % legend on first plot only
        leg = legend([labels string2latex(yM_labels(i)) string2latex(y_labels(i))], ...
            'Interpreter', 'Latex', 'Location', leg_loc);
    end
    %set(leg, 'Position', [0.415 0.762 0.2624 0.1449])
    txt = sprintf('(%s) Observer estimates of process output', char(96+i));
    if ny > 1
        txt = strjoin({txt, char(string(i))});
    end
    title(txt, 'Interpreter','Latex')
    grid on
end

idx = find(~u_known);
for i = 1:nw
    axs(ny+i) = subplot(ny+nw,1,ny+i);
    labels = {};
    for j = obs_sel
        plot(t_range2, X_est(:, idx(i) + n*(j-1)), 'Linewidth', 2); hold on
        label = sprintf("$%s$ by %s", x_est_plot_labels{idx(i)}, ...
            escape_latex_chars(obs_labels{j}));
        labels = [labels {label}];
    end
    stairs(t_range2, Pd(:, i), 'k-'); hold on;
    xlim(t_range2([1 end])')
    ylim(y_lims{i_plot}); i_plot = i_plot + 1;
    set(gca, 'TickLabelInterpreter', 'latex')
    if i == nw  % xlabel on last plot only
        xlabel(time_label, 'Interpreter', 'Latex');
    end
    y_label = string2latex(strjoin([p_labels(i) x_est_plot_labels(idx(i))], ', '));
    ylabel(y_label,'Interpreter','Latex')
    if i == 1  % legend on first plot only
        leg = legend([labels string2latex(p_labels(i))], ...
            'Interpreter', 'Latex', 'Location', 'best');
    end
    %set(leg, 'Position', [0.415 0.29 0.2624 0.1449])
    txt = sprintf('(%s) Observer estimates of input disturbance %d', char(96+i+ny));
    if nw > 1
        txt = strjoin({txt, char(string(i))});
    end
    title(txt, 'Interpreter','Latex')
    grid on
end

linkaxes(axs, 'x')

% Standard size is [560  420]
set(gcf, 'Position', [50 130 448 168*(ny+nw)])

filename = sprintf('%s_y_est2', sim_name);
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Plot cumulative output estimation errors over time

% Select time range
t_stop = 5000;  % e.g. 300 or 5000

% Get output estimation errors for all observers
t_range = sim_results{t <= t_stop, 't'};
E_obs = sim_results{t <= t_stop, ...
    startsWith(sim_results.Properties.VariableNames, 'E_obs')};

% Resample while keeping peak values
nT_min = 500;
[~, ~, Esq_max_obs] = downsample_evenly(E_obs(1:end-1,:).^2, nT_min);

% Calculate cumulative sum of squared errors
E_cum_obs = cumsum(E_obs.^2, 1);

% Resample using mean values
[E_cum_obs, ~, ~] = downsample_evenly(E_cum_obs(1:end-1,:), nT_min);

%nT = size(E_cum_obs, 1);
%if nT > 600
%    n_step = floor(nT / 500);
%end
%t_range = t_range(1:n_step:end);
% Resample using min values
[~, t_range, ~] = downsample_evenly(t_range(1:end-1), nT_min);

switch ny

    case 1

        figure(30); clf

        % Plot squared errors
        ax1 = subplot(2,1,1);
        make_tsplot(Esq_max_obs(:, (obs_sel-1)*ny+i), t_range, ...
            escape_latex_chars(obs_labels(obs_sel)), time_label)
        if i == 1
            xlabel([])
        end
        y_err_label = strjoin({y_labels{i}, y_est_plot_labels{i}}, '-');
        y_sqerr_label = strjoin(["(" y_err_label ")^2"], "");
        ylabel(string2latex(y_sqerr_label), 'Interpreter','Latex')
        title('(a) Squared output estimation errors', 'Interpreter', 'Latex')

        % Plot cumulative squared errors
        ax2 = subplot(2,1,2);
        make_tsplot(E_cum_obs(:, (obs_sel-1)*ny+i), t_range, ...
            escape_latex_chars(obs_labels(obs_sel)), time_label)
        y_label = string2latex(strjoin(["\sum_1^k{" y_sqerr_label "}"], ""));
        ylabel(string2latex(y_label), 'Interpreter', 'Latex')
        title('(b) Cumulative sum of squared estimation errors', 'Interpreter', 'Latex')
        legend(ax1, 'off');
    
        linkaxes([ax1, ax2], 'x')

        % Standard size is [560  420]
        set(gcf, 'Position', [610 950-(410*(i-1)) 448 336])
    
        filename = sprintf('%s_cum_err_y%d', sim_name, i);
        save_fig_to_pdf(fullfile(plot_dir, filename))

    otherwise

        figure(30); clf

        axs = zeros(1, ny);
        for i = 1:ny

            axs(i) = subplot(ny,1,i);
            make_tsplot(E_cum_obs(:, (obs_sel-1)*ny+i), t_range, ...
                escape_latex_chars(obs_labels(obs_sel)), time_label)
            y_err_label = strjoin({y_labels{i}, y_est_plot_labels{i}}, '-');
            y_label = string2latex(strjoin({'\sum_1^k{(', y_err_label, ')^2}'}, ''));
            ylabel(string2latex(y_label), 'Interpreter', 'Latex')
            txt = sprintf('(%s) Cumulative sum of squared errors - output %d', char(96+i), i);
            title(txt, 'Interpreter', 'Latex')
            if i > 1
                legend(axs(i), 'off');
            end
            if i < ny
                xlabel("");
            end    
        end
        linkaxes(axs, 'x')
    
        % Standard size is [560  420]
        set(gcf, 'Position', [610 950-(410*(i-1)) 448 336])
    
        filename = sprintf('%s_cum_err_y%d', sim_name, i);
        save_fig_to_pdf(fullfile(plot_dir, filename))

end


%% Bar plot summary of RMSEs of state estimates

figure(5); clf

for i = 1:n

    X_rmse_i_labels = compose('RMSE_x_est_%d_%s', i, string(obs_sel_labels));

    subplot(n,1,i)
    y1_values = results_summary_no_noise{results_summary_no_noise.seed == seed, X_rmse_i_labels};
    y2_values = results_summary{results_summary.seed == seed, X_rmse_i_labels};
    plot_labels = ["no noise", "due to noise"];
    if any(y1_values > y2_values)
        warning("RMSE values without noise higher")
        assert(all(y1_values/y2_values < 1.0001)) % sometimes KF1 is slightly better with noise
        y1_values(y1_values > y2_values) = y2_values(y1_values > y2_values);
    end

    x_labels = categorical(escape_latex_chars(obs_sel_labels));
    x_labels = reordercats(x_labels, string(x_labels));
    bar(x_labels, [y1_values; y2_values-y1_values], 'stacked');
    set(gca, 'TickLabelInterpreter', 'latex')
    ylim(axes_limits_with_margin(y2_values, 0.25, [0 0.01]))
    if i == n
        xlabel('Observer', 'Interpreter', 'Latex')
    end
    y_axis_label = "$\mathrm{RMSE}(\hat{X}_" + string(i) ...
        + ",X_" + string(i) + ")$";
    ylabel(y_axis_label, 'Interpreter', 'Latex')
    grid on
    if i == 1
        legend(plot_labels, 'Interpreter', 'Latex')
    end
    title_txt = sprintf("(%s) Root-mean-squared errors of state %d estimates", char(96+i), i);
    title(title_txt, 'Interpreter', 'Latex')

    % Display RMSEs
    fprintf("RMSEs of state %d\n", i)
    disp(array2table([y1_values; y2_values-y1_values], ...
    'VariableNames', obs_sel_labels, 'RowNames', {'Without noise', 'With noise'}))

end

% Standard size is [560  420]
set(gcf, 'Position', [1170 400 420 168*n])

filename = sprintf('%s_x_err_bar', sim_name);
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Bar plot summary of RMSEs of output estimates

% Make labels
switch ny
    case 1
        Y_rmse_i_labels = compose('RMSE_y_est_%s', string(obs_sel_labels));
        y_axis_labels = "$\mathrm{RMSE}(\hat{Y},Y)$";
    case 2
        Y_rmse_i_labels = cell(ny, numel(obs_sel_labels));
        y_axis_labels = cell(1, ny);
        for i = 1:ny
            Y_rmse_i_labels(i, :) = compose('RMSE_y_est_%d_%s', i, string(obs_sel_labels));
            y_axis_labels{i} = sprintf('$%s{RMSE}(%s{Y}_%d,Y_%d)$', '\mathrm', '\hat', i, i);
        end
        y_axis_labels = string(y_axis_labels);
end

figure(6); clf

for i = 1:ny

    subplot(ny,1,i)
    y1_values = results_summary_no_noise{results_summary_no_noise.seed == seed, Y_rmse_i_labels(i, :)};
    y2_values = results_summary{results_summary.seed == seed, Y_rmse_i_labels(i, :)};
    if any(y1_values > y2_values)
        warning("RMSE values without noise higher")
        assert(all(y1_values/y2_values < 1.0001)) % sometimes KF1 is slightly better with noise
        y1_values(y1_values > y2_values) = y2_values(y1_values > y2_values);
    end
    x_labels = categorical(escape_latex_chars(obs_sel_labels));
    x_labels = reordercats(x_labels, string(x_labels));
    bar(x_labels, [y1_values; y2_values-y1_values], 'stacked');
    set(gca, 'TickLabelInterpreter', 'latex')
    ylim(axes_limits_with_margin(y2_values, 0.25, [0 0.01]))
    if i == ny
        xlabel('Observer', 'Interpreter', 'Latex')
    end
    ylabel(y_axis_labels(i), 'Interpreter', 'Latex')
    grid on
    if i == 1
        legend("no noise", "due to noise", 'Interpreter', 'Latex')
    end
    if ny > 1
        title_txt = sprintf("(%s) Root-mean-squared errors of output %d estimates", char(96+i), i);
        title(title_txt, 'Interpreter', 'Latex')
    end

    % Display RMSEs
    if ny > 1
        fprintf("RMSEs of output %d\n", i);
    else
        fprintf("RMSEs of output\n");
    end
    disp(array2table([y1_values; y2_values], ...
    'VariableNames', obs_sel_labels, 'RowNames', {'Without noise', 'With noise'}))

end

% Standard size is [560  420]
set(gcf, 'Position', [1170 900 420 168*ny])

filename = sprintf('%s_y_err_bar', sim_name);
save_fig_to_pdf(fullfile(plot_dir, filename))


%% Box plot of RMSE results

% Plots distribution of results over all seeds and all output
% estimates (if multiple simulations exist).

% Note: exclude the results for the simulation outputs used to
% select the observer parameters

seed_values = unique(results_summary.seed);
n_seeds = numel(seed_values);
if n_seeds > 4

    assert(n_seeds == numel(results_summary.seed))
    if ny == 1
        seed_sel = ~(seed_values == 6);
    elseif ny == 2
        seed_sel = ~(seed_values == 0);
    end
    seed_sel_values = seed_values(seed_sel);
    fprintf("Seeds to include: %s\n", strjoin(string(seed_sel_values)));

    RMSE_data = nan(sum(seed_sel), n_obs, ny);
    for i = 1:ny
        RMSE_data(:, :, i) = results_summary{seed_sel, Y_rmse_i_labels(i, :)};
    end
    RMSE_means = mean(RMSE_data, 3);

    figure(7); clf

    boxplot(RMSE_means, 'Labels', escape_latex_chars(obs_sel_labels))
    xlabel('Observer', 'Interpreter', 'Latex')
    if ny == 1
        y_axis_label = "$\mathrm{RMSE}(\hat{Y},Y)$";
    else
        y_axis_label = "$\mathrm{RMSE}(\hat{\mathbf{Y}},\mathbf{Y})$";
    end
    ylabel(y_axis_label, 'Interpreter', 'Latex')
    set(gca, 'TickLabelInterpreter', 'latex')
    grid on

    % Standard size is [560  420]
    set(gcf, 'Position', [1170 100 420 180])
    
    filename = sprintf('%s_y_err_box', sim_name);
    save_fig_to_pdf(fullfile(plot_dir, filename))

end

% Display summary table
mean_RMSE_means = mean(RMSE_means);
mean_RMSE_mean_KF3 = mean_RMSE_means(find(obs_sel_labels == "KF3"));
rmse_summary = array2table( ...
    [mean_RMSE_means; ...
     mean_RMSE_means./mean_RMSE_mean_KF3-1], ...
    'VariableNames', obs_sel_labels, ...
    'RowNames', {'mean RMSE', 'RMSE / RMSE KF3 -1'} ...
);
disp(rmse_summary)


return


%% Plot MKF filter results
%
% These were copied from file 
%  - robertson-1995/rod_obs_sim_analysis_MKF.m
% But it doesn't work yet.


% Choose MKF observer to make detailed plots for
obs_label = "MKF_SF95";


%% Waterfall plot of trace of P

trP = sim_results_obs{:, compose("trP_%d", 1:n_filt)};

figure(10); clf
ax_labels = {'Time ($t$)', 'Filter $f$', '$\mathrm{Tr}(P_f(k))$'};
make_waterfall_plot(t(t_sel), trP(t_sel, 1:min(n_filt, n_max)), [0 2.2], ...
    ax_labels, view_angles);

set(gcf, 'Position', [610, 925, 450, plot_height])
filename = sprintf("%s_%s_trP_wfplot", sim_name, obs_label);
save_fig_to_pdf(fullfile(plot_dir, filename));


%% Waterfall plot of conditional probabilities

p_seq_g_Yk = sim_results_obs{:, compose("p_seq_g_Yk_%d", 1:n_filt)};

figure(11); clf
ax_labels = {'Time ($t$)', 'Filter $f$', '$Pr(\Gamma_f(k) \mid Y(k))$'};
make_waterfall_plot(t(t_sel), p_seq_g_Yk(t_sel, 1:min(n_filt, n_max)), [0 1], ...
    ax_labels, view_angles);

set(gcf, 'Position', [610, 600, 450, plot_height])
filename = sprintf("%s_%s_p_seq_g_Yk_wfplot", sim_name, obs_label);
save_fig_to_pdf(fullfile(plot_dir, filename));


%% Plot MKF observer filter estimates

% Number of model states
n = sim_summary{1, 'n'};

% True state values
X = sim_results{:, compose("X_%d", 1:n)};

% MKF observer estimates
X_est = sim_results{:, compose("X_est_%d", (i_obs-1)*n+1:i_obs*n)};

% MKF observer filter estimates
X_est_f = sim_results_obs{:, compose("X_est_%d", 1:n_filt*n)};

figure(12); clf

make_MKF_X_est_plot(X_est(t_sel, :), t(t_sel), X_est_f(t_sel, :), X(t_sel, :))

set(gcf, 'Position', [610, 175, 450, 336])
filename = sprintf("%s_%s_x_est", sim_name, obs_label);
save_fig_to_pdf(fullfile(plot_dir, filename));


%% Plot observer shock sequences (Gamma(k))
% obs_seq = repmat(cell2mat(obs.seq)', ceil(nT+1 / (obs.d * obs.f)), 1);
% obs_seq = obs_seq(1:nT+1, :);
% 
% figure(12); clf
% p1 = show_waterfall_plot(t(151:251), obs_seq(151:251, :), [0 1], ...
%     ax_labels, [0 82]);
% hold on
% p2 = show_waterfall_plot(t(151:251), p_seq_g_Yk(151:251, :), [0 1], ...
%     ax_labels, [0 82], filepath);

