function make_MKF_X_est_plot(X_est,t,X_est_f,X,nh_max)
% Makes a figure showing time-series subplots of the
% estimates produced by a mult-model observer (MKF).
% Each subplot shows the overall observer estimate,
% the estimates by each of its filters, and (optionally)
% the true state values.

    if nargin < 5
        % Maximum number of filters to list in legend
        nh_max = 5;
    end
    plot_X = nargin >= 4;
    n = size(X_est, 2);
    assert(rem(size(X_est_f, 2), n) == 0)
    n_filt = size(X_est_f, 2) / n;
    nf_show = min(n_filt, nh_max);
    nT = size(X_est, 1);
    if nT <= 100
        plot_func = @stairs;
    else
        plot_func = @plot;
    end

    % Loop over each state variable and make subplot
    for i = 1:n

        y_label = sprintf("$%s_{%d}(k|k), %s_{%d,f}(k|k), x_{%d}(k)$", ...
            "\hat{x}", i, "\hat{x}", i, i);
    
        subplot(n,1,i)
    
        plot_func(t, X_est(:, i), 'Linewidth', 2); hold on
        plot_func(t, X_est_f(:, i:n:n_filt*n))
        if plot_X
            plot_func(t, X(:, i), 'k--')
        end
        ylim(axes_limits_with_margin([X_est(:, i) X_est_f(:, i:n:n_filt*n)]));
        grid on
        if i == n
            xlabel('Time ($t$)')
        end
        ylabel(y_label)
        title_text = sprintf("(%s) Estimates of model state $x_{%d}(k)$", char(96+i), i);
        title(title_text)

        % Make array of filter labels for legend
        si = sprintf("%d", i);
        x_est_f_labels = "$\hat" + ...
            compose(strcat("{x}_{", si, ",%d}"), 1:nf_show) + "$";
        labels = [{"$\hat"+sprintf('{x}_{%d}$', i)} x_est_f_labels];
        if plot_X
            labels = [labels {sprintf('$x_{%d}$', i)}];
        end
        legend(labels, 'Interpreter', 'latex', 'Location', 'best');

    end

end