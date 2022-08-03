function make_tsstatplot(Y, t, y_labels, x_label, y_lim, area, line)
% make_tsstatplot(Y, t, y_labels, x_label, y_lim, area, line)
% Time-series plots of the statistics of groups of signals.
%
% Arguments
%   Y : array of group data or cell array of group data
%     arrays.
%   t : column vector of time values.
%   y_labels : label or cell array of labels for each data
%     group.
%   x_label : x-axis label (optional, default is '$t$')
%   y_lim : y-axis limits (optional, default is nan(2))
%   area : char indicating the statistics represented by
%     the fill area. Options include:
%      - 'minmax' for minimum and maximum values (this is
%        the default),
%      - 'std' for -/+ one standard deviation from mean,
%      - string representation of a percentile, e.g. '90'
%        for the 5 to 95 percent range.
%   line : char indicating the statistic represented by
%     the solid line. Options include:
%      - 'mean' for the mean value (this is the default),
%      - 'median' for the median value.
%
    if nargin < 7
        line = 'mean';
    end
    if nargin < 6
        area = 'minmax';
    end
    if nargin < 5
        y_lim = nan(2);
    end
    if nargin < 4
        x_label = "$t$";
    else
        x_label = string(x_label);
    end
    if isnumeric(Y)  % case of only one data group
        Y = {Y};
    end
    if nargin < 3
        if numel(Y) == 1
            y_labels = "$y(t)$";
        else
            y_labels = compose("$y_{%d}(t)$", 1:numel(Y));
        end
    else
        y_labels = string(y_labels);
    end
    line_labels = cell(1, numel(y_labels)*2);
    % Get color order
    colors = get(gca,'colororder');
    set(gca, 'ColorOrder', colors);
    for iy = 1:numel(Y)
        switch area
            case 'minmax'
                Y_upper = max(Y{iy}, [], 2);
                Y_lower = min(Y{iy}, [], 2);
                fill_labels = 'min, max';
            case 'std'
                Y_avg = mean(Y{iy}, 2);
                Y_std = std(Y{iy}, [], 2);
                Y_upper = Y_avg + Y_std;
                Y_lower = Y_avg - Y_std;
                fill_labels = '+/- 1 std. dev.';
            otherwise
                pct = str2num(area);
                if isequaln(pct, [])
                    error("ValueError: invalid line type")
                end
                assert((pct > 0) & (pct < 100))
                pcts = [(100-pct)*0.5 100-(100-pct)*0.5];
                Y_lower_upper = prctile(Y{iy}', pcts)';
                Y_upper = Y_lower_upper(:, 2);
                Y_lower = Y_lower_upper(:, 1);
                fill_labels = char(sprintf('$%d\\%s$ CI', pct, '%'));
        end
        switch line
            case 'mean'
                Y_line = mean(Y{iy}, 2);
            case 'median'
                Y_line = median(Y{iy}, 2);
            otherwise
                error("ValueError: invalid line type")
        end
        line_labels{iy*2-1} = [y_labels{iy} ' ' fill_labels];
        line_labels{iy*2} = [y_labels{iy} ' ' line];
        % Modify colors if plotting more than one group
        if numel(Y) > 1
            colors = get(gca,'colororder');
        end
        % Make filled area plot
        inBetween = [Y_lower; flip(Y_upper)];
        t2 = [t; flip(t)];
        fill(t2, inBetween, colors(iy, :), 'LineStyle', 'none'); 
        alpha(.25); hold on
        % Add line plot
        h = plot(t, Y_line, 'Linewidth', 2);
        %set(h, {'color'}, {colors(1, :); colors(2, :)});
        set(h, {'color'}, {colors(iy, :)});
    end
    ylim(axes_limits_with_margin([Y_upper Y_lower], 0.1, y_lim, y_lim))
    set(gca, 'TickLabelInterpreter', 'latex')
    if strlength(x_label) > 0
        xlabel(x_label, 'Interpreter', 'Latex')
    end
    ylabel(strjoin(y_labels, ', '), 'Interpreter', 'latex')
    legend(line_labels, 'Interpreter', 'latex', 'Location', 'best')
    grid on
