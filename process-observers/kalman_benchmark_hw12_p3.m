% Kalman filter simulation code from GEL 7029 Course, Homework 12.
% See file /gel-7029/homework/hw12/hw12_p3_kalman.m

clear all; close all; clc;

% Create directory for plot output
if ~exist('plots','dir')
    mkdir('plots')
end

% Initialize random number generator
rng(0);

% Load system from file
sys_test_siso2

% Observer parameters
W = 0.5; % estimate of Wp used in the filter design
V = 0.8; % estimate of Vp used in the filter design
P = eye(n)*1000; % Initialize covariance matrix
Q = diag(repmat(W,n,1));
R = diag(repmat(V,ny,1));
N = zeros(n,ny);

% Steady-state Kalman filter
Gkf = eye(n);
Hkf = zeros(ny,n);
Gmodel = ss(A,[B Gkf],C,[D Hkf],Ts);
[KF_ss, K_ss, P_ss] = kalman(Gmodel,Q,R,N,'delayed');

% Simulation settings
nT = 100; % number of timesteps
xp = [0; 0]; % initial state vector of the plant
xe = [0; 0]; % initial estimates

% seed random number generator
rng(0)

% Input signal - pseudo-random binary sequence
warnId = 'Controllib:estimation:initialPRBSSequence';
warnStruct = warning('off',warnId);
u = idinput(nT,'PRBS',[0 0.5]);
warning(warnStruct);

% Measurement noise for the whole simulation
v = sqrt(Rp)*randn(nT,1);

% Process noise for the whole simulation
w = (sqrt(Qp)*randn(2,nT))';

% Run simulation with inputs u, v, w
[t, x_p, x_est, y_p, y_est, k_est, p_est] = ...
    run_simulation(u,v,w,A,B,C,P,Q,R,xp,xe,Ts);

% Estimation errors
e_est = x_p - x_est;
% Get diagonal elements of covariance matrix
p_diag = [p_est(:,1,1) p_est(:,2,2)];

sim_results = table( ...
    t(1:nT), u, v, w, x_p, x_est, e_est, y_p, y_est, k_est, p_diag, ...
    'variableNames', {'t', 'u(t)', 'v(t)', 'w(t)', 'X(t)', 'X_e(t)', 'E(t)', ...
                      'y(t)', 'y_est(t)', 'K(t)', 'P(t)'} ...
);

% Display first 10 rows of simulation results
head(sim_results, 10)

% Compare to steady-state Kalman filter
assert(max(abs(k_est(end-1,:)' - K_ss)) < 1e-6)
assert(max(abs(p_est(end-1,:)' - P_ss(:))) < 1e-6)

% Save simulation results for unit testing
filename = 'hw12_p3_kalman_sim_benchmark.csv';
writetable(sim_results, fullfile('results', filename))


% %% Figures showing the time evolution of different signals
% % Plot comparing true and estimated states
% figure
% for i=1:n
%     subplot(n,1,i)
%     plot(t(1:nT),x_p(:,i),'.-'); hold on
%     plot(t(1:nT),x_est(:,i),'.-')
%     xlabel('t')
%     title(sprintf('Plant state and observer estimate - x_%d',i))
%     legend('State', 'Estimate')
%     grid on
%     set(gca,'FontSize',12);
% end
% set(gcf,'position',[100,100,600,400]);
% saveas(gcf,'plots/hw12_3_states.png')

return


%% Kalman filter update equations
function [K, P] = kalman_update(P,A,C,Q,R)
    % Compute observer gains
    K = A*P*C'*(C*P*C' + R)^-1;
    % Update covariance matrix
    P = A*(P - P*C'*(C*P*C' + R)^-1*C*P)*A' + Q;
end


%% Simulation function
function [t, x_p, x_est, y_p, y_est, k_est, p_est] = ...
         run_simulation(u,v,w,A,B,C,P,Q,R,xp,xe,Ts)
    nT = size(u,1);
    n = length(A); % number of states
    ny = size(C,1); % number of outputs
    % Matrices to record the results
    t = [0:nT]'*Ts;
    x_p = zeros(nT, n); % Plant states
    y_p = zeros(nT, ny); % Plant outputs
    x_est = zeros(nT, n); % State estimates
    y_est = zeros(nT, ny); % Output estimates
    k_est = zeros(nT, n); % Observer gains
    p_est = zeros(nT, n, n); % Covariance matrix diagonals

    % Simulation loop
    for i = 1:nT
        % Compute plant output
        y = C*xp + v(i); % plant output
        % Save current states and estimates
        x_p(i,:) = xp;
        x_est(i,:) = xe;
        y_p(i,:) = y;
        y_est(i,:) = C*xe;
        % Compute observer gains and update covariance matrix
        [K, P] = kalman_update(P,A,C,Q,R);
        % Compute next state estimates
        xe = A*xe + B*u(i) + K*(y-C*xe);
        % Save data
        p_est(i,:,:) = P;
        k_est(i,:) = K';
        % Update process states for next iteration
        xp = A*xp + B*u(i) + w(i,:)'; % plant state update
    end
    % Save final states and estimates
    x_p(i,:) = xp;
    x_est(i,:) = xe;
    y_p(i,:) = y;
    y_est(i,:) = C*xe;
end
