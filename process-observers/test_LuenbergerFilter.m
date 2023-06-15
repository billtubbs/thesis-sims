% Test LuenbergerFilter class
% (See class definition in LuenbergerFilter.m)

clear all

% Simulation settings
nT = 200; % number of points
Ts = 3; % sampling period
t = Ts*(0:nT)'; % time vector

% SISO system example from GEL-7029 in file Luenb_no_integr.mlx
A = [0.82 0; 0 0.9];
B = [1; -1];
C = [-0.5 1];
D = 0;
Qp = diag([0.3; 0.2]); % covariance - process noise
Rp = 0.4; % variance - measurement noise

% Dimensions
n = size(A, 1); % number of states
nu = size(B, 2);  % number of inputs
ny = size(C, 1); % number of outputs

% Check if benchmark simulation data file exists
if ~isfile('results/Luenberger_Filter_sim_benchmark.csv')
    error("Run 'Luenb_no_integr_benchmark.mlx' to generate benchmark data.")
end


%% Define and simulate steady-state Kalman filters

% Covariance matrices
Q = diag([0.1; 0.1]);
R = 0.5;
N = zeros(n,ny);

% Prepare a struct of model parameters
model = struct();
model.A = A;
model.B = B;
model.C = C;
model.Ts = Ts;
model.Q = Q;
model.R = R;

% Define Luenberger observer
label = "LB1";
poles = [0.9; 0.9];
x0 = [0.1; 0.5];
LB = LuenbergerFilter(model,poles,label,x0);
assert(strcmp(LB.type, "LB"))
assert(isequal(LB.xkp1_est, x0))
assert(LB.ykp1_est == C * x0)
assert(LB.label == label);
assert(LB.n == n)
assert(LB.nu == nu)
assert(LB.ny == ny)

% Re-define with no initial state specified (should be set to zero)
LB = LuenbergerFilter(model,poles,label);
assert(isequal(LB.xkp1_est, zeros(n, 1)))
assert(LB.ykp1_est == 0)
K_test = [0.16; 0];
assert(isequal(round(LB.K, 4), K_test))
assert(isequal(LB.xkp1_est, zeros(2, 1)))
assert(LB.ykp1_est == 0)

% seed random number generator
rng(0)

% Measurement noise for the whole simulation
v = sqrt(Rp)*randn(nT,1);

% Process noise for the whole simulation
w = sqrt(Qp)*randn(2,nT);

% Output disturbance
p = zeros(nT,1);
p(t>=300) = 1; % step at time t=300

u0 = 1;  % initial value of u
x0 = inv(eye(length(A)) - A) * B * u0;  % steady-state value of x

% Intialize system (at k = 0)
x = x0;

% Input signal
U = [zeros(10,1); ones(nT+1-10, 1)]; % process input for the whole simulation

% Matrices to collect simulation data
xNprocess = zeros(n, nT+1); % process states
yNprocess = zeros(ny, nT+1); % process outputs
xNobserver = zeros(n, nT+1); % estimated states
yNobserver = zeros(ny, nT+1); % estimated process outputs
xNkalman2 = zeros(n, nT+1); % estimated states
yNkalman2 = zeros(ny, nT+1); % estimated process outputs

for i = 1:nT

    % Process output in current timestep
    y = C*x + v(i) + p(i);
    
    % Record process states and output
    xNprocess(:, i) = x;
    yNprocess(:, i) = y; 

    % Process states in next timestep
    x = A*x + B*U(i) + w(:,i);

    % Lunberger filter update
    LB.update(y, U(i));

    % Record Kalman filter estimates in next timestep
    xNobserver(:, i+1) = LB.xkp1_est;
    yNobserver(:, i+1) = LB.ykp1_est;

end
t = Ts * (0:nT)';


% plot results

% figure(1); clf
% 
% subplot(411);
% plot(t', yNprocess, 'k', t', yNkalman1, 'r', t', yNkalman2, 'g', 'Linewidth', 2)
% legend('Process output', 'KF1 estimates', 'KF2 estimates')
% ylabel('y_1')
% grid on
% 
% subplot(412);
% plot(t', xNprocess(1,:), 'k', t', xNkalman1(1,:), 'r', ...
%     t', xNkalman2(1,:), 'g', 'Linewidth', 2)
% legend('Process state', 'KF1 estimates', 'KF2 estimates')
% ylabel('x_1')
% grid on
% 
% subplot(413);
% plot(t', xNprocess(2,:), 'k', t', xNkalman1(2,:), 'r', ...
%     t', xNkalman2(2,:), 'g', 'Linewidth', 2);
% legend('Process state', 'KF1 estimates', 'KF2 estimates')
% ylabel('x_2')
% grid on
% 
% subplot(414);
% stairs(t', U', 'Linewidth', 2);
% xlabel('Time [s]');
% ylabel('u_1')
% grid on

% Display results
sim_results = [table(t,U) ...
    array2table(xNprocess', 'VariableNames', {'x1', 'x2'}) ...
    array2table(xNobserver', 'VariableNames', {'x1_est', 'x2_est'}) ...
    array2table(yNprocess', 'VariableNames', {'y'}) ...
    array2table(yNobserver', 'VariableNames', {'y_est'}) ...
];

%head(sim_results)

% Verify results by comparing with Luenb_no_integr_benchmark.mlx

filename = 'Luenberger_Filter_sim_benchmark.csv';
bench_sim_results = readtable(fullfile('results', filename));

%head(bench_sim_results)

% Check states
assert(isequal( ...
    round(sim_results{1:200, {'x1', 'x2'}}, 6), ...
    round(bench_sim_results{1:200, {'x1', 'x2'}}, 6) ...
))

% Check state estimates
assert(isequal( ...
    round(sim_results{1:200, {'x1_est', 'x2_est'}}, 6), ...
    round(bench_sim_results{1:200, {'x1_est', 'x2_est'}}, 6) ...
))
