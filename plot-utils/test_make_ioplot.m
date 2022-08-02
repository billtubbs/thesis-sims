% Test function make_ioplot.m

clear all; close all

% Directory to save test plots
plot_dir = 'plots';
if ~isfolder(plot_dir)
    mkdir(plot_dir)
end


%% IO plot of discrete-time SISO system

Ts = 0.5;
nT = 10;
t = Ts*(0:nT)';
U = zeros(nT+1,1);
U(t >= 1, :) = 1;
G = tf(1, [1 1], Ts);
[Y, t] = lsim(G,U,t);

% Basic plot with default labels
figure(1)
make_ioplot(Y, t, U)
save_fig_to_pdf(fullfile(plot_dir, 'ioplot1_def.pdf'))

% With custom labels
u_labels = {'$u(t)$'};
y_labels = {'$y(t)$'};
x_label = '$t$ (seconds)';
titles_text = {'(a) Outputs', '(b) Inputs'};

figure(2)
make_ioplot(Y, t, U, u_labels, y_labels, x_label, nan(2), nan(2), ...
    titles_text, 'stairs')
save_fig_to_pdf(fullfile(plot_dir, 'ioplot1.pdf'))

% Same plot with strings and string arrays
u_labels = "$u(t)$";
y_labels = "$y(t)$";
titles_text = ["(a) Outputs" "(b) Inputs"];
x_label = "$t$ (seconds)";

figure(3)
make_ioplot(Y, t, U, u_labels, y_labels, x_label, nan(2), nan(2), ...
    titles_text, 'stairs')
save_fig_to_pdf(fullfile(plot_dir, "ioplot1s.pdf"))

% With axes limits
figure(4)
x_label = '$t$ (seconds)';
y1_lim = [-1 1];
y2_lim = [-1 1];
make_ioplot(Y, t, U, u_labels, y_labels, x_label, y1_lim, y2_lim, ...
    {'', ''}, 'stairs')
save_fig_to_pdf(fullfile(plot_dir, 'ioplot2.pdf'))


%% IO plot of continuous-time 2x2 system

% Simulate continuous-time 2x2 system
t = linspace(0, 10, 101)';
nT = size(t, 1) - 1;
U = zeros(nT+1,2);
U(t >= 1, 1) = 1;
U(t >= 3, 2) = -1;
G = [tf(1, [1 1]) 0;
     0            tf(1, [2 1])];
[Y, t] = lsim(G,U,t);

% Basic plot with default labels
figure(5)
make_ioplot(Y, t, U)
save_fig_to_pdf(fullfile(plot_dir, 'ioplot3_def.pdf'))

% With custom labels
u_labels = {'$u_1(t)$', '$u_2(t)$'};
y_labels = {'$y_1(t)$', '$y_2(t)$'};

figure(6)
make_ioplot(Y, t, U, u_labels, y_labels)
save_fig_to_pdf(fullfile(plot_dir, 'ioplot3.pdf'))

% Same plot with strings and string arrays
u_labels = ["$u_1(t)$" "$u_2(t)$"];
y_labels = ["$y_1(t)$" "$y_2(t)$"];

figure(7)
make_ioplot(Y, t, U, u_labels, y_labels)
save_fig_to_pdf(fullfile(plot_dir, 'ioplot3s.pdf'))

figure(8)
x_label = '$t$ (seconds)';
y1_lim = [-1.5 1.5];
y2_lim = [-1.5 1.5];
make_ioplot(Y, t, U, u_labels, y_labels, x_label, y1_lim, y2_lim)
save_fig_to_pdf(fullfile(plot_dir, 'ioplot4.pdf'))


%% Discretre-time 4x4 system

Ts = 1;
nT = 50;
t = Ts*(0:nT)';
U = zeros(nT+1, 4);
U((t >= 2 & t < 25), 1) = 1;
U(t >= 5, 2) = 1;
U(t >= 10, 3) = -1;
U(t >= 15, 4) = -1;

G = [tf(1, [1 1]) tf(-0.2, [1 1]) 0 0;
     tf(-0.2, [1 1]) tf(1, [1 1]) 0 0;
     0 0 tf(-1, [1 1]) tf(0.5, [1 1]);
     0 0 tf(0.5, [1 1]) tf(-1, [1 1])];
Gd = c2d(G, Ts);
[Y, t] = lsim(G,U,t);

% Basic plot with default labels
figure(9)
make_ioplot(Y, t, U)
save_fig_to_pdf(fullfile(plot_dir, 'ioplot5_def.pdf'))

% With custom labels
u_labels = {'$u_1(t)$', '$u_2(t)$', '$u_3(t)$', '$u_4(t)$'};
y_labels = {'$y_1(t)$', '$y_2(t)$', '$y_2(t)$', '$y_4(t)$'};

figure(10)
make_ioplot(Y, t, U, u_labels, y_labels)
save_fig_to_pdf(fullfile(plot_dir, 'ioplot5.pdf'))

close all