function make_iodmplot(Y, Y_m, Y_model, t, U, u_labels, y_labels, ...
    x_label, y1_lim, y2_lim, titles_text, kind, intr)
% make_iodmplot(Y, Y_m, Y_model, t, U, u_labels, y_labels, ...
%     x_label, y1_lim, y2_lim, titles_text, kind, intr)
% Time series plots of input signal, true output signal, measured
% output, and model predicted output data.

    if nargin < 13
        intr = 'latex';
    end
    if nargin < 12
        kind = 'plot';
    end
    if nargin < 11
        titles_text = ["(a) Outputs" "(b) Inputs"];
    else
        titles_text = string(titles_text);
    end
    if nargin < 10
        y2_lim = nan(2);
    end
    if nargin < 9
        y1_lim = nan(2);
    end
    if nargin < 8
        x_label = "$t$";
    else
        x_label = string(x_label);
    end
    if nargin < 7
        if size(Y, 2) == 1
            y_labels = "$y(t)$";
        else
            y_labels = compose("$y_{%d}(t)$", 1:size(Y, 2));
        end
    else
        y_labels = string(y_labels);
    end
    if nargin < 6
        if size(U, 2) == 1
            u_labels = "$u(t)$";
        else
            u_labels = compose("$u_{%d}(t)$", 1:size(U, 2));
        end
    else
        u_labels = string(u_labels);
    end
    switch kind
        case 'plot'
            plotf = @plot;
        case 'stairs'
            plotf = @stairs;
    end

    ax1 = subplot(3,1,1:2);
    ny = size(Y_m, 2);
    colors = get(ax1, 'ColorOrder');
    if ny > 1
        colors(ny+1:ny*2, :) = colors(1:ny, :);
    end
    if size(Y_m, 1) < 100
        d_style = 'o';
    else
        d_style = '.';
    end
    set(ax1, 'ColorOrder', colors);
    plotf(t, Y, 'Linewidth', 2); hold on
    set(ax1, 'ColorOrder', colors);
    plotf(t, Y_m, d_style)
    set(ax1, 'ColorOrder', colors);
    plotf(t, Y_model, 'Linewidth', 2)
    ylim(axes_limits_with_margin([Y Y_m], 0.1, y1_lim, y1_lim))
    set(gca, 'TickLabelInterpreter', intr)
    if numel(y_labels) > 3
        y_axis_labels = [y_labels(1), y_labels(2), {'...'}];
        ylabel(strjoin(y_axis_labels, ', '), 'Interpreter', intr)
    else
        ylabel(strjoin(y_labels, ', '), 'Interpreter', intr)
    end
    if numel(y_labels) > 1
        legend(y_labels, 'Interpreter', intr, 'Location', 'best')
    end
    if titles_text{1}
        title(titles_text{1}, 'Interpreter', intr)
    end
    grid on

    ax2 = subplot(3,1,3);
    make_tsplot(U, t, u_labels, x_label, y2_lim, titles_text{2}, 'stairs')

    linkaxes([ax1, ax2], 'x')