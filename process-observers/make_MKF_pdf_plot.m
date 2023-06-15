function make_MKF_pdf_plot(obs, f, yk, y_lim, p_max)
% Display probability density of the output estimate
% of filter f of the multi-model observer compared to
% the current data point.
%

    if nargin < 4
        y_lim = [-inf inf];
    end

    % Get y_est(k/k-1) estimated in previous time step
    %TODO: This will involve storing yk_pred for each observer
    yk_est = obs.filters{f}.ykp1_est;

    % Calculate covariance of the output estimation errors
    P = obs.filters{f}.P;
    C = obs.filters{f}.C;
    yk_cov = C*P*C' + obs.filters{f}.R;

    % Make sure covariance matrix is symmetric
    if ~isscalar(yk_cov)
        yk_cov = triu(yk_cov.',1) + tril(yk_cov);
    end

    x_label = "$y(k),\hat{y}(k),y_M(k)$";
    labels = {'$p(y(k))$', '$\hat{y}(k)$', '$y_m(k)$'};
    
    % Alternatively, to calc a y range based on std. dev.s
    % from mean use:
    % sd = sqrt(yk_cov); y_lim = y_mean + 3*[-sd sd];

    make_pdf_plot(yk_est, yk_cov, yk, labels, x_label, y_lim, p_max)

end