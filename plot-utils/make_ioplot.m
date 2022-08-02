function make_ioplot(Y, t, U, u_labels, y_labels, x_label, y1_lim, ...
    y2_lim, titles_text, kind, intr)
% make_ioplot(Y, t, U, u_labels, y_labels, x_label, y1_lim, ...
%     y2_lim, titles_text, kind)
% Time series plots of input and output signals

    if nargin < 11
        intr = 'latex';
    end
    if nargin < 10
        kind = 'plot';
    end
    if nargin < 9
        titles_text = ["(a) Outputs" "(b) Inputs"];
    else
        titles_text = string(titles_text);
    end
    if nargin < 8
        y2_lim = nan(2);
    end
    if nargin < 7
        y1_lim = nan(2);
    end
    if nargin < 6
        x_label = "$t$";
    else
        x_label = string(x_label);
    end
    if nargin < 5
        if size(Y, 2) == 1
            y_labels = "$y(t)$";
        else
            y_labels = compose("$y_{%d}(t)$", 1:size(Y, 2));
        end
    else
        y_labels = string(y_labels);
    end
    if nargin < 4
        if size(U, 2) == 1
            u_labels = "$u(t)$";
        else
            u_labels = compose("$u_{%d}(t)$", 1:size(U, 2));
        end
    else
        u_labels = string(u_labels);
    end

    ax1 = subplot(2,1,1);
    make_tsplot(Y, t, y_labels, [], y1_lim, titles_text(1), kind, intr)

    ax2 = subplot(2,1,2);
    make_tsplot(U, t, u_labels, x_label, y2_lim, titles_text(2), ...
        'stairs', intr)

    linkaxes([ax1, ax2], 'x')
