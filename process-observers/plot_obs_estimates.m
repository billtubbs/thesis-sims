function plot_obs_estimates(t,X,X_est,Y,Y_est,obs_labels,intp)
% Display plots of observer estimates compared to
% true values.
    if nargin < 7
        intp = 'latex';
    end
    n = size(X, 2);
    ny = size(Y, 2);
    axs = nan(1, ny+n);

    for i = 1:ny
        y_values = [Y(:,i) Y_est(:, i:ny:end)];
        axs(i) = subplot(ny+n,1,i);
        plot(t, y_values(:, 2:end),'Linewidth',2); hold on
        plot(t, y_values(:, 1), 'k--');
        max_min = [min(y_values, [], [1 2]) max(y_values, [], [1 2])];
        bd = max([0.1 diff(max_min)*0.1]);
        ylim(max_min + [-bd bd])
        y_label = sprintf("$y_%d(k)$", i);
        ylabel(y_label, 'Interpreter', intp)
        title(strjoin(["Output" y_label]),'Interpreter',intp)
        if nargin > 5
            obs_labels = string(obs_labels);
            legend([obs_labels "true"], 'Interpreter', intp, ...
                'Location', 'best')
        end
        grid on
    end

    for i = 1:n
        y_values = [X(:,i) X_est(:, i:n:end)];
        axs(ny+i) = subplot(ny+n,1,ny+i);
        plot(t, y_values(:, 2:end), 'Linewidth', 2); hold on
        plot(t, y_values(:, 1),'k--');
        max_min = [min(y_values, [], [1 2]) max(y_values, [], [1 2])];
        bd = max([0.1 diff(max_min)*0.1]);
        ylim(max_min + [-bd bd])
        y_label = sprintf("$x_%d(k)$", i);
        ylabel(y_label,'Interpreter',intp)
        if nargin > 5
            legend([obs_labels "true"], 'Interpreter', intp, ...
                'Location', 'best')
        end
        title(strjoin(['State' y_label]), 'Interpreter', intp)
        grid on
    end
    xlabel('Time', 'Interpreter', intp)

    linkaxes(axs, 'x')

end