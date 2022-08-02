% Test function make_iorplot.m

clear all; close all

% Directory to save test plots
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir)
end


%% discrete-time SISO system

Ts = 0.5;
nT = 10;
t = Ts*(0:nT)';
U = zeros(nT+1,1);
U(t >= 1, :) = 1;
G = tf(1, [1 1], Ts);
[Y, t] = lsim(G,U,t);
R = ones(nT+1, 1);
u_labels = {'$u(t)$'};
y_labels = {'$y(t)$'};
r_labels = {'$r(t)$'};

figure(1)
titles_text = {'(a) Outputs', '(b) Inputs'};
make_iorplot(Y, t, U, R, u_labels, y_labels, r_labels, '$t$', ...
    nan(2), nan(2), titles_text, 'stairs')
save_fig_to_pdf(fullfile(plot_dir, 'iorplot1.pdf'))

figure(2)
x_label = '$t$ (seconds)';
y1_lim = [-1 1];
y2_lim = [-1 1];
make_iorplot(Y, t, U, R, u_labels, y_labels, r_labels, x_label, ...
    y1_lim, y2_lim, {'', ''}, 'stairs')
save_fig_to_pdf(fullfile(plot_dir, 'iorplot1.pdf'))


%% Continuous-time 2x2 system
t = linspace(0, 10, 101)';
nT = size(t, 1) - 1;
U = zeros(nT+1,2);
U(t >= 1, 1) = 1;
U(t >= 3, 2) = -1;
G = [tf(1, [1 1]) 0;
     0            tf(1, [2 1])];
[Y, t] = lsim(G,U,t);
R = zeros(nT+1, 2);
R(t >= 1, 1) = 1;
R(t >= 3, 2) = -1;
u_labels = string2latex({'u_1(t)', 'u_2(t)'});
y_labels = string2latex({'y_1(t)', 'y_2(t)'});
r_labels = string2latex({'r_1(t)', 'r_2(t)'});

figure(3)
make_iorplot(Y, t, U, R, u_labels, y_labels, r_labels)
save_fig_to_pdf(fullfile(plot_dir, 'iorplot3.pdf'))
