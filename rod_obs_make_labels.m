%% Generates all the labels used in tables and plots
%
% Inputs:
%  - system parameters: n, ny, nu, nw
%  - observers : cell array containing observers

% Make signal labels
u_labels = vector_element_labels('u', '(k)', nu, false);
x_labels = vector_element_labels('x', '(k)', n, false);
p_labels = vector_element_labels('p', '(k)', nw, false);
y_labels = vector_element_labels('y', '(k)', ny, false);
y_m_labels = vector_element_labels('y_m', '(k)', ny, false);
y_model_labels = vector_element_labels('y_{model}', '(k)', ny, false);

% Combine all input signal labels
input_labels = repmat({''}, 1, nu);
input_labels(u_meas') = u_labels;
input_labels(~u_meas') = p_labels;

% Make observer labels into an array
n_obs = numel(observers);
obs_labels = @(i) observers{i}.label;
obs_labels = cellfun(obs_labels, num2cell(1:n_obs), ...
    'UniformOutput', false);

% Make array of x_est(k), y_est(k) labels
x_est_labels = vector_element_labels('x_est', '', n, false);
x_est_plot_labels = vector_element_labels('\hat{x}', '(k|k-1)', n, false);
y_est_labels = vector_element_labels('y_est', '', ny, false);
y_est_plot_labels = vector_element_labels('\hat{y}', '(k|k-1)', ny, false);
y_m_plot_labels = vector_element_labels('y_m', '(k)', ny, false);

% Make metrics labels for all observers, e.g. for observer 'KF1':
%  - 'MSE_y_est_KF1' : overall MSE
%  - 'MSE_tr_y_est_KF1' : MSE in transition periods
%  - 'MSE_ss_y_est_KF1' :  MSE in steady-state periods
%  - 'Var_ss_y_est_KF1' : Variance in steady-state periods
metrics_labels = {'MSE', 'MSE_tr', 'MSE_ss', 'Var_ss', 'MSD_ss'};
n_metrics = numel(metrics_labels);
obs_metrics_labels = cell(n_metrics, n_obs * ny);
for i = 1:n_metrics
    metric_label = metrics_labels{i};
    labels = matrix_element_labels(metric_label, y_est_labels, obs_labels, '');
    obs_metrics_labels(i, :) = labels(:)';
end
