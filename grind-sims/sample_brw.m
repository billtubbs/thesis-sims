%% Bounded Random Walk Disturbance Model
%
% Author: Bill Tubbs, April 2022.
%
% Simulation of the Bounded Random Walk stochastic process proposed
% by J. Nicolau:
% Title: Stationary Processes That Look Like Random Walks - The Bounded
%        Random Walk Process in Discrete and Continuous Time
% Author: J. Nicolau
% Publication: Econometric Theory, 18, 2002, 99-118.
%

clear all

rng(0)
plot_dir = "plots";

% Model parameters
sd_e = 0.4;    % Noise std. dev.
beta = -15;   % k parameter in Nicolau's paper
alpha1 = 3;
alpha2 = 3;
tau = 100;   % initial value x(0)

% Simulation timesteps
N = 1000;

% Simulation
p = sample_bounded_random_walk(sd_e, beta, alpha1, alpha2, N, tau);

figure(1); clf
t = 0:N-1;
plot(t, p, 'Linewidth', 2);
ylim(axes_limits_with_margin([p; 95; 100], 0.05));
yline(95, 'k--');
yline(105, 'k--');
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$p(k)$', 'Interpreter', 'latex')
grid on
legend('BRW', 'Interpreter', 'latex', 'location', 'best')
set(gcf,'Position',[100 100 400 150])
save2pdf(fullfile(plot_dir, 'brw_sim.pdf'))


%% Plot of Bounded Random Walk bias function

figure(2); clf
x = 94:0.1:106;
% a = brw_reversion_bias(x, 0, alpha2, beta, tau);
% plot(x, a, 'Linewidth', 2); hold on
% a = brw_reversion_bias(x, alpha1, 0, beta, tau);
% plot(x, a, 'Linewidth', 2);
a = brw_reversion_bias(x, alpha1, alpha2, beta, tau);
plot(x, a, 'Linewidth', 2);
set(gca, 'TickLabelInterpreter', 'latex')
grid on
xlabel("$p(k)$", 'Interpreter', 'latex')
ylabel("$a(p(k))$", 'Interpreter', 'latex')
%title("BRW Reversion Function")
set(gcf,'Position',[100 300 300 150])
save2pdf(fullfile(plot_dir, 'brw_a.pdf'))


%% Stationary pdf

dt = 0.02;
x = 94:dt:106;
brw_pdf = brwpdf(x, alpha1, alpha2, beta, tau, sd_e);

% Normalize the pdf
brw_pdf = brw_pdf ./ (sum(brw_pdf) * dt);

figure(3); clf
plot(x, brw_pdf, 'Linewidth', 2)
set(gca, 'TickLabelInterpreter', 'latex')
grid on
xlabel("$p(k)$", 'Interpreter', 'latex')
ylabel("$\mathrm{Pr}(p(k))$", 'Interpreter', 'latex')
ylim(axes_limits_with_margin(brw_pdf));
set(gcf,'Position',[100 525 300 150])
%save2pdf(fullfile(plot_dir, 'brw_pdf.pdf'))


%% Comparison with RW pdf

N = 20;
rw_pdf = normpdf(x, tau, N(1)*sd_e^2);

figure(6); clf
plot(x, rw_pdf, 'Linewidth', 2); hold on
plot(x, brw_pdf, 'Linewidth', 2)
set(gca, 'TickLabelInterpreter', 'latex')
grid on
xlabel("$p(k)$", 'Interpreter', 'latex')
ylabel("$\mathrm{Pr}(p(k))$", 'Interpreter', 'latex')
ylim(axes_limits_with_margin([rw_pdf' brw_pdf']));
labels = {sprintf('RW ($k = %d$)', N), 'BRW ($k = \infty$)'};
legend(labels, 'Interpreter', 'latex', 'location', 'south')
set(gcf,'Position',[100 750 300 150])
save2pdf(fullfile(plot_dir, 'brw_pdf.pdf'))
