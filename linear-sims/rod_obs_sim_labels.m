% Prepares label strings for data storage and plotting

addpath("~/ml-plot-utils/")

% Make signal labels
y_labels = vector_element_labels('y', '(k)', ny, false);
switch ny
    case 1
        yM_labels = {'y_M(k)'};
    otherwise
        yM_labels = compose('y_{M,%d}(k)', 1:ny);
end
u_labels = vector_element_labels('u', '(k)', sum(u_known), false);
x_labels = vector_element_labels('x', '(k)', n, false);
p_labels = vector_element_labels('p', '(k)', nw, false);

% Combine all input signal labels
input_labels = repmat({''}, 1, nu);
input_labels(u_known') = u_labels;
input_labels(~u_known') = p_labels;

% Make array of x_est(k), y_est(k) labels
x_est_labels = vector_element_labels('x_est', '', n, false);
x_est_plot_labels = vector_element_labels('\hat{x}', '(k|k)', n, false);
y_est_labels = vector_element_labels('y_est', '', ny, false);
y_est_plot_labels = vector_element_labels('\hat{y}', '(k|k)', ny, false);

% Make array of RMSE labels for each observer
RMSE_x_est_labels = matrix_element_labels('RMSE', x_est_labels, obs_labels, '');
RMSE_x_est_labels = RMSE_x_est_labels(:)';
RMSE_y_est_labels = matrix_element_labels('RMSE', y_est_labels, obs_labels, '');
RMSE_y_est_labels = RMSE_y_est_labels(:)';

% X-axis label for time-series plots
time_label = 'Time ($t$)';