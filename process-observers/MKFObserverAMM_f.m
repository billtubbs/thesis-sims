% Multi-model Kalman Filter class definition
%
% obs = MKFObserverAMM_f(A,B,C,Ts,P0,Q,R,label,x0,p_seq_g_Yk_init)
% Class for simulating an autonomous multi-model (AMM)
% Kalman filter for state estimation of system with more than
% one possible mode (multiple models). This is used as the base
% class for the following multiple model observers:
%
%  - MKFObserverGPB1 - 1st order Generalised Pseudo-Bayes algorithm
%  - MKFObserverGPB2 - 2nd order Generalised Pseudo-Bayes algorithm
%  - MKFObserverIMM - Interacting multiple model algorithm.
% 
% These observers are in 'filtering form', which means they
% produce posterior estimates of the states and outputs at the 
% current time given the data at the  current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next time 
% instant given the data at the current time are also calculated:
%
%  xkp1_hat(k+1|k) : estimate of states at time k + 1
%  ykp1_hat(k+1|k) : estimate of outputs at time k + 1
%
% The system model is defined as:
%
%   x(k+1) = A(k) x(k) + B(k) u(k) + w(k)
%     y(k) = C(k) x(k) + v(k)
%
% Note: It is assumed that there is no direct transmission
% (D = 0).
%
% Arguments:
%   A, B, C : Cell arrays containing discrete-time system
%       matrices for each switching system modelled.
%   Ts : Sample period.
%   P0 : Initial covariance matrix of the state estimates
%       (same for each filter).
%   Q : Cell array of process noise covariance matrices for
%       each switching system.
%   R : Cell array of output measurement noise covariance
%       matrices for each switching system.
%   label : string name.
%   x0 : Initial state estimates (optional, default zeros).
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%

classdef MKFObserverAMM_f < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        Ts (1, 1) double {mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
        n_filt (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        A cell
        B cell
        C cell
        D cell
        Q cell
        R cell
        label (1, 1) string
        P0 double
        x0 (:, 1) double
        gamma_k double {mustBeInteger, mustBeNonnegative}
        p_seq_g_Yk_init double
        p_seq_g_Ykm1 double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_gammak_g_Ykm1 double
        p_gamma_k_g_gamma_km1 double
        filters  cell
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        Pkp1 double
        ykp1_est (:, 1) double
        type (1, 1) string
    end
    methods
        function obj = MKFObserverAMM(A,B,C,Ts,P0,Q,R,label,x0, ...
                p_seq_g_Yk_init)

            % System dimensions
            [n, nu, ny] = check_dimensions(A{1}, B{1}, C{1});

            % Number of switching systems
            nj = numel(A);

            % Number of filters required (= no. of models in
            % this case)
            n_filt = nj;

            if nargin < 10
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(n_filt, 1) ./ double(n_filt);
            end
            if nargin < 9
                x0 = zeros(n,1);
            end
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.Ts = Ts;
            obj.Q = Q;
            obj.R = R;
            obj.label = label;
            obj.P0 = P0;
            obj.x0 = x0;

            % Model indicator values (zero-based)
            obj.gamma_k = (0:nj-1)';
            obj.p_seq_g_Yk_init = p_seq_g_Yk_init;

            % Check all other system matrix dimensions have same 
            % input/output dimensions and number of states.
            for j = 2:nj
                [n_j, nu_j, ny_j] = check_dimensions(A{j}, B{j}, C{j});
                assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
                    "ValueError: size of A, B, and C")
            end

            % Store other useful variables
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;
            obj.nj = nj;
            obj.n_filt = n_filt;
            obj.type = "MKF_AMM";

            % Initialize all variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Initial values of prior conditional probabilities at k = -1 
            obj.p_seq_g_Yk = obj.p_seq_g_Yk_init;

            % Empty vectors to store values for filter calculations
            % p(y(k)|Gamma(k),Y(k-1))
            obj.p_yk_g_seq_Ykm1 = nan(obj.n_filt, 1);
            % Pr(gamma(k)|Y(k-1))
            obj.p_gammak_g_Ykm1 = nan(obj.n_filt, 1);
            % Pr(Gamma(k)|Gamma(k-1)) = 1 since no transitions with AMM
            obj.p_gamma_k_g_gamma_km1 = ones(obj.nj, 1);
            % Pr(Gamma(k)|Y(k-1))
            obj.p_seq_g_Ykm1 = nan(obj.n_filt, 1);

            % Instantiate Kalman filters for each model
            obj.filters = cell(obj.n_filt, 1);
            fmt = strcat('%s%0', ...
                char(string(strlength(string(obj.n_filt)))), 'd');
            % Initialize each filter
            for i = 1:obj.n_filt
                label_i = sprintf(fmt,obj.label,i);
                % Index of system model
                ind = obj.gamma_k(i) + 1;
                obj.filters{i} = KalmanFilterF(obj.A{ind},obj.B{ind}, ...
                    obj.C{ind},obj.Ts,obj.P0, obj.Q{ind},obj.R{ind}, ...
                    label_i,obj.x0);
            end
            %TODO: Could just reset the KFs using their reset methods.

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.C{1} * obj.xkp1_est;

            % At initialization at time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);

            % Initialize error covariance at k = 0
            obj.Pk = nan(obj.n);
            obj.Pkp1 = obj.P0;

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
        % updates the multi-model Kalman filter and calculates the
        % estimates of the states and output at the next sample
        % time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %

            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

            % Arrays to collect estimates from each filter
            Xkf_est = zeros(obj.n, 1, obj.n_filt);
            Pkf_est = zeros(obj.n, obj.n, obj.n_filt);
            Ykf_est = zeros(obj.ny, 1, obj.n_filt);
            Xkp1f_est = zeros(obj.n, 1, obj.n_filt);
            Pkp1f_est = zeros(obj.n, obj.n, obj.n_filt);
            Ykp1f_est = zeros(obj.ny, 1, obj.n_filt);

            % Bayesian update to conditional probabilities
            for f = 1:obj.n_filt

                % Update observer estimates, gain and covariance matrix
                obj.filters{f}.update(yk, uk);
                assert(~any(isnan(obj.filters{f}.xkp1_est)))

                % Compute posterior probability density of y(k)
                % using the posterior PDF (assumed to be a normal 
                % distribution) and output estimate.

                % Get y_f_est(k/k-1) estimated in previous time step
                ykf_est = obj.filters{f}.ykp1_est;

                % Covariance of the output estimation errors
                Sk = obj.filters{f}.Sk;

                % Make sure covariance matrix is symmetric
                if ~isscalar(Sk)
                    Sk = triu(Sk.',1) + tril(Sk);
                end

                % Calculate normal probability density (multivariate)
                obj.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, ykf_est, Sk);

                % Save state and output estimates
                Xkf_est(:, :, f) = obj.filters{f}.xk_est';
                Pkf_est(:, :, f) = obj.filters{f}.Pk;
                Ykf_est(:, :, f) = obj.filters{f}.yk_est';
                Xkp1f_est(:, :, f) = obj.filters{f}.xkp1_est';
                Pkp1f_est(:, :, f) = obj.filters{f}.Pkp1;
                Ykp1f_est(:, :, f) = obj.filters{f}.ykp1_est';

            end

            assert(~any(isnan(obj.p_yk_g_seq_Ykm1)))
            assert(~all(obj.p_yk_g_seq_Ykm1 == 0))

            % Compute Pr(Gamma(k)|Y(k-1)) in current timestep from
            % previous estimate (Pr(Gamma(k-1)|Y(k-1))) and transition
            % probabilities
            obj.p_seq_g_Ykm1 = obj.p_gamma_k_g_gamma_km1 .* obj.p_seq_g_Yk;

            % Bayesian update of Pr(Gamma(k)|Y(k))
            likelihood = obj.p_yk_g_seq_Ykm1 .* obj.p_seq_g_Ykm1;
            obj.p_seq_g_Yk = likelihood ./ sum(likelihood);
            % To prevent likelihoods going to zero use this:
            %obj.p_seq_g_Yk = likelihood * 0.998 ./ sum(likelihood) + 0.001;
            % Note: above calculation normalizes p_seq_g_Yk so that
            % assert(abs(sum(obj.p_seq_g_Yk) - 1) < 1e-15) % is always true

            % Compute multi-model observer state and output estimates
            % and estimated state error covariance using the weighted-
            % averages based on the conditional probabilities.
            [obj.xk_est, obj.yk_est, obj.Pk] = ...
                weighted_avg_estimates(Xkf_est, Ykf_est, Pkf_est, ...
                    obj.p_seq_g_Yk);

            % Weighted average estimates at next time instant
            [obj.xkp1_est, obj.ykp1_est, obj.Pkp1] = ...
                weighted_avg_estimates(Xkp1f_est, Ykp1f_est, Pkp1f_est, ...
                    obj.p_seq_g_Yk);

        end
    end
end
