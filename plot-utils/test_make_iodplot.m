% Test functions make_iodplot.m and make_iodmplot.m

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
Y_m = Y + 0.1*randn(size(Y));  % add measurement noise
u_labels = {'$u(t)$'};
y_labels = {'$y(t)$', '$y_m(t)$'};
id_data = iddata(Y_m,U,Ts);
model = polyest(id_data, [0 1 0 0 1 1]);
[Y_model, t] = lsim(model,U,t);

figure(1); clf
titles_text = {'(a) Outputs', '(b) Inputs'};
make_iodplot(Y, Y_m, t, U, u_labels, y_labels, '$t$', nan(2), nan(2), ...
    titles_text, 'stairs')
save_fig_to_pdf(fullfile(plot_dir, 'iodplot1.pdf'))

figure(2); clf
titles_text = {'(a) Outputs', '(b) Inputs'};
y_labels = {'$y(t)$', '$y_m(t)$', '$y_{model}(t)$'};
make_iodmplot(Y, Y_m, Y_model, t, U, u_labels, y_labels, '$t$', ...
    nan(2), nan(2), titles_text, 'stairs')
save_fig_to_pdf(fullfile(plot_dir, 'iodplot2.pdf'))


%% Continuous-time 2x2 system
t = linspace(0, 10, 101)';
nT = size(t, 1) - 1;
U = zeros(nT+1,2);
U(t >= 1, 1) = 1;
U(t >= 3, 2) = -1;
G = [tf(1, [1 1]) 0;
     0            tf(1, [2 1])];
[Y, t] = lsim(G,U,t);
Y_m = Y + 0.1*randn(size(Y));  % add measurement noise
u_labels = string2latex({'u_1(t)', 'u_2(t)'});
y_labels = string2latex({'y_1(t)', 'y_2(t)', 'y_{m,1}(t)', 'y_{m,2}(t)'});

Ts = diff(t(1:2));
id_data = iddata(Y_m,U,Ts);
poly_order = [zeros(2, 2) ones(2, 2) zeros(2, 1) zeros(2, 1) ones(2, 2) ones(2, 2)];
model = polyest(id_data, poly_order);
[Y_model, t] = lsim(model,U,t);

figure(3); clf
x_label = '$t$ (seconds)';
make_iodplot(Y, Y_m, t, U, u_labels, y_labels, x_label)
save_fig_to_pdf(fullfile(plot_dir, 'iodplot3.pdf'))

figure(4); clf
y_labels = string2latex({'y_1(t)', 'y_2(t)', 'y_{m,1}(t)', 'y_{m,2}(t)', ...
    'y_{model,1}(t)', 'y_{model,2}(t)'});
make_iodmplot(Y, Y_m, Y_model, t, U, u_labels, y_labels, '$t$ (seconds)')
save_fig_to_pdf(fullfile(plot_dir, 'iodplot4.pdf'))
