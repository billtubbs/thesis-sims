% Test GPB update functions derived from code by Zhao et al.
%
% Function:
%  - GPB2_estimation.m
%
% These functions were adapted from code provided in:
%  - S. Zhao, C. K. Ahn, P. Shi, Y. S. Shmaliy, and F. Liu, (2019).
%     Bayesian State Estimations for Markovian Jump Systems.
%

clear
rng(0)

% Use this if making plots
%addpath("../plot-utils")

data_dir = 'data';
if ~isfolder(data_dir)
    mkdir(data_dir)
end

% Define switching system
nj = 2;  % number of modes
A = cell(1, nj);
B = cell(1, nj);
C = cell(1, nj);
D = cell(1, nj);

% Sampling period
Ts = 0.5;

% Discrete time state space models
% Model #1
A{1} = 0.7;
B{1} = 1;
C{1} = 0.3;
D{1} = 1;  % Note: D(r(k)) is part of noise model, not direct transmission

% Model #2
A{2} = 0.9;
B{2} = 1;
C{2} = -0.3;  % -ve gain!
D{2} = 1;

% Transition probabilities
T = [0.95 0.05; 0.01 0.99];
assert(all(sum(T, 2) == 1))

% Dimensions
n = size(A{1}, 1);
nu = size(B{1}, 2);
ny = size(C{1}, 1);

% Struct to store variables for GPB1 algorithm
GPB1 = struct();

% Define observer model
GPB1.Model.TP = T;
GPB1.Model.A = A;
GPB1.Model.Bu = B;
GPB1.Model.Bw = {eye(n), eye(n)};
GPB1.Model.C = C;
GPB1.Model.D = D;

% Struct to store variables for GPB2 algorithm
GPB2 = GPB1;

% Observer parameters
P0 = 0.1;
Q1 = 0.1;
Q2 = Q1;  % same for both modes
R1 = 0.1^2;
R2 = R1;  % same for both modes
GPB1.Model.Q = {Q1, Q2};
GPB1.Model.R = {R1, R2};
GPB2.Model.Q = {Q1, Q2};
GPB2.Model.R = {R1, R2};

% Initial conditions
x0 = 0;

% Initialise GPB1 variables
GPB1.x = x0;
GPB1.P = P0;
GPB1.mu = [0.5; 0.5];

% Initialise GPB2 variables
GPB2.x = repmat(x0, 1, nj);
GPB2.P = repmat(P0, 1, 1, nj);
GPB2.mu = [0.5; 0.5];

% Measurement noise std. dev.
sigma_M = 0.1;  % TODO: greater than 0.2 causes numerical errors

% Simulation settings
nT = 60;
t = Ts*(0:nT)';

% Inputs
U = zeros(nT+1,1);
U(t > 2) = 1;
V = sigma_M * randn(size(t));

% Actual system mode sequence
seq = int8(ones(nT+1, 1));
seq(t >= 10, 1) = 2;

% Arrays to store simulation results
X = nan(nT+1, n);
Y = nan(nT+1, ny);
Ym = nan(nT+1, ny);
Xk_est.GPB1 = nan(nT+1, n);
Yk_est.GPB1 = nan(nT+1, ny);
P_modes.GPB1 = nan(nT+1, nj);
Xk_est.GPB2 = Xk_est.GPB1;
Yk_est.GPB2 = Yk_est.GPB1;
P_modes.GPB2 = P_modes.GPB1;

% Run simulation
xk = x0;
for i = 1:nT+1

    % Inputs
    uk = U(i);
    vk = V(i);

    % Calculate system states and output in current timestep
    rk = seq(i);
    xk = A{rk} * xk + B{rk} * uk;
    yk = C{rk} * xk;
    ymk = yk + D{rk} * vk;

    % Store system state and outputs
    X(i, :) = xk';
    Y(i, :) = yk';
    Ym(i, :) = ymk';

    % Update GPB1 estimator
    [GPB1.x,GPB1.P,GPB1.mu] = ...
        GPB1_estimation(GPB1.x,GPB1.P,ymk,uk,GPB1.Model,GPB1.mu);

    assert(~any(isnan(GPB1.P), [1 2]))

    % Store state estimates
    Xk_est.GPB1(i, :) = GPB1.x;
    P_modes.GPB1(i, :) = GPB1.mu';

    % Calculate y_est(k)
    Yk_est.GPB1(i, :) = sum( ...
        [GPB1.Model.C{1}*GPB1.x; GPB1.Model.C{2}*GPB1.x] .* GPB1.mu ...
    )';

    % Update GPB2 estimator
    [GPB2.x,GPB2.P,GPB2.mu,Out_x,Out_P] = ...
        GPB2_estimation(GPB2.x,GPB2.P,ymk,uk,GPB2.Model,GPB2.mu);

    assert(~any(isnan(GPB2.P), [1 2 3]))

    % Store state estimates
    Xk_est.GPB2(i, :) = Out_x';
    P_modes.GPB2(i, :) = GPB2.mu';

    % Calculate y_est(k)
    Yk_est.GPB2(i, :) = sum( ...
        [GPB2.Model.C{1}*GPB2.x(:, 1); GPB2.Model.C{2}*GPB2.x(:, 2)] ...
            .* GPB2.mu ...
    )';

end

% Calculate output estimation errors
rmses.GPB1 = sqrt(mean((Yk_est.GPB1 - Y).^2));
rmses.GPB2 = sqrt(mean((Yk_est.GPB2 - Y).^2));

switch sigma_M
    case 0.1
        assert(round(rmses.GPB1, 6) == 0.070806)
        assert(round(rmses.GPB2, 6) == 0.073353)
    case 0.2
        assert(round(rmses.GPB1, 6) == 0.143784)
        assert(round(rmses.GPB2, 6) == 0.154601)
    otherwise
        warning("RMSE not checked")
end

% Save results to file
sim_results = table(t,U,seq,X,Y,V,Ym,Yk_est.GPB1,Yk_est.GPB2,Xk_est.GPB1,Xk_est.GPB2);
sim_results.Properties.VariableNames = { ...
    't', 'U', 'seq', 'X', 'Y', 'V', 'Ym', ...
    'Yk_est_GPB1','Yk_est_GPB2','Xk_est_GPB1','Xk_est_GPB2' ...
};
filename = "GPB_sim_results.csv";
writetable(sim_results,fullfile(data_dir,filename))
% This data is used to test MKFObserverGPB1 and MKFObserverGPB2


% % Plot of system inputs and outputs
% figure(1); clf
% make_iodplot(Y,Ym,t,[U seq],{'$u(k)$','$r(k)$'},{'$y(k)$','$y_m(k)$'})
% 
% % Plot of observer estimates compared to true outputs
% figure(2); clf
% plot_obs_estimates(t,X,[Xk_est.GPB1 Xk_est.GPB2],Y, ...
%     [Yk_est.GPB1 Yk_est.GPB2],{'GPB1','GPB2'})
% 
% % Plot of estimated mode probabilities GPB1
% figure(3); clf
% subplot(2,1,1)
% bar(t, P_modes.GPB1, 'stacked', 'BarWidth', 1);
% ylim([0 1])
% set(gca, 'TickLabelInterpreter','latex')
% xlabel("Time ($t$)",'Interpreter','latex')
% ylabel("$Pr(r(k))$",'Interpreter','latex')
% legend("$Pr(r(k)=1)$","$Pr(r(k)=2)$",'Interpreter','latex')
% title("Estimated mode probabilities - GPB1",'Interpreter','latex')
% 
% subplot(2,1,2)
% bar(t, [seq==1 seq==2], 'stacked', 'BarWidth', 1);
% ylim([0 1])
% set(gca, 'TickLabelInterpreter','latex')
% xlabel("Time ($t$)",'Interpreter','latex')
% ylabel("$Pr(r(k))$",'Interpreter','latex')
% legend("$r(k)=1$","$r(k)=2$",'Interpreter','latex')
% title("True system mode",'Interpreter','latex')
% 
% % Plot of estimated mode probabilities GPB2
% figure(4); clf
% subplot(2,1,1)
% bar(t, P_modes.GPB2, 'stacked', 'BarWidth', 1);
% ylim([0 1])
% set(gca, 'TickLabelInterpreter','latex')
% xlabel("Time ($t$)",'Interpreter','latex')
% ylabel("$Pr(r(k))$",'Interpreter','latex')
% legend("$Pr(r(k)=1)$","$Pr(r(k)=2)$",'Interpreter','latex')
% title("Estimated mode probabilities - GPB2",'Interpreter','latex')
% 
% subplot(2,1,2)
% bar(t, [seq==1 seq==2], 'stacked', 'BarWidth', 1);
% ylim([0 1])
% set(gca, 'TickLabelInterpreter','latex')
% xlabel("Time ($t$)",'Interpreter','latex')
% ylabel("$Pr(r(k))$",'Interpreter','latex')
% legend("$r(k)=1$","$r(k)=2$",'Interpreter','latex')
% title("True system mode",'Interpreter','latex')