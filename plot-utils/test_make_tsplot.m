% Test function make_tsplot.m

clear all; close all

% Directory to save test plots
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir)
end


%% Time series plot of single discrete signal

Ts = 0.5;
nT = 10;
t = Ts*(0:nT)';
X = randn(nT+1, 1);
y_labels = "$Y(t)$";

% Basic plot with default labels
figure(1)
make_tsplot(X, t, y_labels)
save_fig_to_pdf(fullfile(plot_dir, 'tsplot1_def.pdf'))

% With custom x-axis labels
x_label = "$t$ (seconds)";
figure(2)
make_tsplot(X, t, y_labels, x_label)
save_fig_to_pdf(fullfile(plot_dir, 'tsplot1a.pdf'))

% With custom title
x_label = '$t$ (seconds)';
title_text = "Time series plot of Y";
figure(3)
make_tsplot(X, t, y_labels, x_label, [nan nan], title_text)
save_fig_to_pdf(fullfile(plot_dir, 'tsplot1b.pdf'))

% With y-axis limits
y_lims = [-1 1];
title_text = 'Time series plot of Y';
figure(4)
make_tsplot(X, t, y_labels, x_label, y_lims, title_text)
save_fig_to_pdf(fullfile(plot_dir, 'tsplot1c.pdf'))

close all


%% Time series plot of multiple signals

Ts = 5;
nT = 10;
t = Ts*(0:nT)';
X = cumsum(randn(nT+1, 4));
y_labels = ["$Y_1(t)$" "$Y_2(t)$" "$Y_3(t)$" "$Y_4(t)$"];

% Basic plot with default labels
figure(5)
make_tsplot(X, t, y_labels)
save_fig_to_pdf(fullfile(plot_dir, 'tsplot2_def.pdf'))

close all
