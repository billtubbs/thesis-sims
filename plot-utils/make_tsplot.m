function make_tsplot(Y, t, y_labels, x_label, y_lim, title_text, kind, ...
    intr)
% make_tsplot(Y, t, y_labels, x_label, y_lim, title_text, kind, intr)
% Time series plots of one or more signals.

    if nargin < 8
        intr = 'latex';
    end
    if nargin < 7
        kind = 'plot';
    end
    if nargin < 6
        title_text = "";
    else
        title_text = string(title_text);
    end
    if nargin < 5
        y_lim = nan(2);
    end
    if nargin < 4
        x_label = "$t$";
    else
        x_label = string(x_label);
    end
    y_labels = string(y_labels);
    switch kind
        case 'plot'
            plotf = @plot;
        case 'stairs'
            plotf = @stairs;
    end
    plotf(t, Y, 'Linewidth', 2);
    ylim(axes_limits_with_margin(Y, 0.1, y_lim, y_lim))
    set(gca, 'TickLabelInterpreter', intr)
    if strlength(x_label) > 0
        xlabel(x_label, 'Interpreter', 'Latex')
    end
    if numel(y_labels) > 3
        y_axis_labels = [y_labels(1), {'...'}, y_labels(end)];
        ylabel(strjoin(y_axis_labels, ', '), 'Interpreter', intr)
    else
        ylabel(strjoin(y_labels, ', '), 'Interpreter', intr)
    end
    if numel(y_labels) > 1
        legend(y_labels, 'Interpreter', intr, 'Location', 'best')
    end
    if strlength(title_text) > 0
        title(title_text, 'Interpreter', intr)
    end
    grid on
