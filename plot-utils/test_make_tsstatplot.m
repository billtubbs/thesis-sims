% Test make_tsstatplot.m

clear all; close all

rng(0);

% Directory to save test plots
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir)
end


%% Example 1 - std. devs. and mean

t = 0.5*(0:20)';
Y = randn(21, 10);

figure(1)
make_tsstatplot(Y, t, 'Output', 'Time (mins)', nan(2), 'std', 'mean')
save_fig_to_pdf(fullfile(plot_dir, 'tsstatplot1.pdf'))

% With y-axis limits
figure(2)
make_tsstatplot(Y, t, '$y(t)$', '$t$ (mins)', [-1 1], 'std', 'mean')
save_fig_to_pdf(fullfile(plot_dir, 'tsstatplot2.pdf'))


%% Example 2 - Confidence intervals and median

% Generate 10 random signals
t = 0.5*(0:20)';
Y = randn(21, 10);

figure(3)
make_tsstatplot(Y, t, '$y(t)$', '$t$ (mins)', nan(2), '90', 'median')
save_fig_to_pdf(fullfile(plot_dir, 'tsstatplot3.pdf'))


%% Example 3 - With more than one group of trajectories

% Generate two groups of 10 random signals
t = 0.5*(0:20)';
Y = {randn(21, 10), randn(21, 10)+2*sin(t)};

figure(4)
y_labels = {'$V(t)$', '$W(t)$'};
make_tsstatplot(Y, t, y_labels, '$t$ (mins)', nan(2), 'minmax', 'mean')
save_fig_to_pdf(fullfile(plot_dir, 'tsstatplot4.pdf'))
