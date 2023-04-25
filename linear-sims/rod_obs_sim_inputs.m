
% Generate input data for observer simulations
%
% Run this script after defining a system using one of
% the following:
%
%  - sys_rodin_step.m
%  - sys_rodin_step_2x2sym.m
%
% Variables which should be defined before running this:
%
%  - A, B, C, D
%  - u_known, y_meas
%  - x0
%  - nT  % number of time steps to simulate
%  - Ts  % sample period
%  - t  % time vector
%  - set the random number seed, e.g. rng(0)
%
% This script is used by run_obs_sim_spec.m during
% automated simulations.
%

addpath("~/process-observers")

% Check necessary variables are already defined
assert(isequal(size(t), [nT+1 1]))
assert(all(diff(t) == Ts))
[n, nu, ny] = check_dimensions(A, B, C, D);
assert(isequal(size(u_known), [nu 1]))
assert(isequal(size(y_meas), [ny 1]))
assert(isequal(size(x0), [n 1]))

% Number of measured and unmeasured input signals
nu = sum(u_known);
nw = sum(~u_known);
assert(isequal(size(p0), [nw 1]))

% Array for combined input signals
U = zeros(nT+1, nu+nw);

switch nu
    case 0  % Arom3 model has no measured inputs

        %TODO: Not implemented yet

    case 1  % SISO system

        % Measured input signals
        U(t >= 5, u_known(1)) = 1;  % add a step

    case 2  % 2x2 system

        % Measured input signals - optional
        %U(t >= 5, u_known(1)) = 1;  % add a step
        %U(t >= 150, u_known(2)) = 1;  % add a step

end

% Number of variances specified for each disturbance input
n_modes = cellfun(@(s) size(s, 2), sigma_wp);

% Only binary modes supported
assert(all(n_modes <= 2))

% Find switching noise models with more than one mode
switching_inputs = find(n_modes > 1);
ns = length(switching_inputs);

% Generate randomly-occurring shock signals
sigma_wp_sw = cell2mat(sigma_wp(switching_inputs)');
[Wp, alpha] = sample_random_shocks([nT+1 ns], epsilon, ...
        sigma_wp_sw(:,2), sigma_wp_sw(:,1));

% Observer model without disturbance noise input
Bu = B(:, u_known);
Du = D(:, u_known);

% Define a multiple-model scheduled Kalman filter based on the
% actual shock signal to use as a benchmark test
switch nw
    case 1

        obs_model = model;  % make copy
        obs_model.B = Bu;
        if isfield(obs_model, "D")
            obs_model.D = Du;
        end
        % P0, R, q1, q2 must all be already defined by sys script
        obs_model.R = R;
        models = {obs_model, obs_model};
        models{1}.Q = diag([q1 sigma_wp_sw(1).^2]);
        models{2}.Q = diag([q1 sigma_wp_sw(2).^2]);
        seq = alpha' + 1;
        SKF = SKFObserverS(models,P0,seq,"SKF");
        observers = [observers {SKF}];

    case 2

        obs_model = model;  % make copy
        obs_model.B = Bu;
        if isfield(obs_model, "D")
            obs_model.D = Du;
        end
        nj = 4;
        % P0, R, q1, q2 must all be already defined by sys script
        obs_model.R = R;
        models = repmat({obs_model}, 1, nj);
        models{1}.Q = diag([q1 q2 sigma_wp_sw(1,1)^2 sigma_wp_sw(2,1)^2]);
        models{2}.Q = diag([q1 q2 sigma_wp_sw(1,2)^2 sigma_wp_sw(2,1)^2]);
        models{3}.Q = diag([q1 q2 sigma_wp_sw(1,1)^2 sigma_wp_sw(2,2)^2]);
        models{4}.Q = diag([q1 q2 sigma_wp_sw(1,2)^2 sigma_wp_sw(2,2)^2]);
        seq = indseq_from_drvs(alpha);
        SKF = SKFObserverS(models,P0,seq,"SKF");
        observers = [observers {SKF}];

    otherwise
        error(">2 RODDs not implemented.")

        % Code from test_MKFObserverSP

%         % Custom MKF test observer
%         % Devise a custom multi-model filter with a shock indicator
%         % sequence that perfectly reflects the shock occurence in
%         % this test simulation (t = t_shock)
%         % Multiple model filter 1
%         A2 = repmat({A}, 1, 3);
%         Bu2 = repmat({Bu}, 1, 3);
%         C2 = repmat({C}, 1, 3);
%         Du2 = repmat({Du}, 1, 3);
%         P0 = 1000*eye(n);
%         %P0_init = repmat({P0}, 1, 3);
%         Q2 = {diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,1)^2]), ...
%               diag([0.01 0.01 sigma_wp(1,2)^2 sigma_wp(2,1)^2]), ...
%               diag([0.01 0.01 sigma_wp(1,1)^2 sigma_wp(2,2)^2])};
%         R2 = {diag(sigma_M.^2), diag(sigma_M.^2), diag(sigma_M.^2)};
%         seq = {zeros(1, nT+1); zeros(1, nT+1); zeros(1, nT+1); zeros(1, nT+1)};
%         seq{2}(t == t_shock(1)) = 1;  % shock 1
%         seq{3}(t == t_shock(2)) = 2;  % shock 2
%         seq{4}(t == t_shock(1)) = 1;  % both
%         seq{4}(t == t_shock(2)) = 2;
%         p_gamma = [1-epsilon epsilon]';
%         Z = [0 0; 1 0; 0 1];  % combinations
%         p_gamma = prod(prob_gamma(Z', p_gamma), 1)';
%         p_gamma = p_gamma ./ sum(p_gamma);  % normalized
%         T = repmat(p_gamma', 3, 1);
%         d = 1;
%         MKF3 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,d,'MKF3');
%         assert(MKF3.n_filt == 4)
%
%         seq = {zeros(1, nT+1)};
%         seq{1}(t == t_shock(1)) = 1;
%         seq{1}(t == t_shock(2)) = 2;
%         MKF4 = MKFObserver(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq,T,d,'MKF4');
%         assert(MKF4.n_filt == 1)
%
%         % Define scheduled Kalman filter
%         % Note: in the case of more than one random input variable, all
%         % possible combinations of the switching systems need to be
%         % accounted for.
%         % Here, we account for 3 possible combinations:
%         % combs = [0 0; 1 0; 0 1];
%         % (This is the same as the MKF filters for the RODD).
%         % seq = sum(alpha .* 2.^(1:-1:0), 2)';
%         SKF = MKFObserverS(A2,Bu2,C2,Du2,Ts,P0,Q2,R2,seq{1},"SKF");

end

% AROM3 model
% % Test with disturbance signals from Fig. 3 of Robertson et al. (1998)
% P = zeros(nT+1, 2);
% step_periods = {t < 4.6, t >= 4.6, t >= 11.2};
% step_values = {[5.4 5.6], [2.7 5.6], [2.7 3]};
% for i = 1:numel(step_periods)
%     P(step_periods{i}, :) = repmat(step_values{i}, sum(step_periods{i}), 1);
% end
% % Shock disturbance that would create this
% Wp = [diff(P - step_values{1}); zeros(1, 2)];
% assert(isequal(size(Wp), [nT+1 nw]))
% % Find when shocks occurred
% [rows,cols,v] = find(Wp);
% alpha = zeros(nT+1, nw);
% for i = 1:numel(rows)
%     alpha(rows(i), cols(i)) = 1;
% end

% % For testing/debugging: Replace above with with this code
% % which sets two input disturbances at fixed times.
% Wp = zeros(nT+1, nw);
% alpha = zeros(size(t));
% t_steps = [100 200];
% Wp(t == t_steps(1), 1) = 1;
% alpha(t == t_steps(1), 1) = 1;
% i = mod(1, nw)+1;
% Wp(t == t_steps(2), i) = -1;
% alpha(t == t_steps(2), i) = 1;

% Combine all model inputs including disturbance signals
U(:, ~u_known) = Wp;

% Measurement noise signals
% TODO: Should measurement noise start at k=0?
V = zeros(nT+1, ny);
V(t >= 0, y_meas) = sigma_M' .* randn(nT+1, ny);

% Process noise signals
W = zeros(nT+1, n);
W(t >= 0, :) = sigma_W' .* randn(nT+1, n);

% Convert disturbance model to state-space
HDdss = ss(HDd);

% Compute disturbance signals at input to plant
% (These are not used by the simulation, only in plots)
[Pd, ~] = lsim(HDdss, Wp, t, p0);

% Combine inputs into a table for use in simulation function
input_data = table(alpha, U, V, W, Pd);
