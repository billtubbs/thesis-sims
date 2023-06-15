% Multi-model Kalman Filter class definition
%
% obs = MKFObserverGPB2(models,P0,T,label,x0,p_seq_g_Yk_init)
%
% Class for simulating the generalised pseudo-Bayes multi-
% model Kalman filter for state estimation of Markov jump
% linear systems. This is the second-order version of the
% algorithm (GPB2).
%
% This is the filtering form of the observer which 
% produces posterior estimates of the states and
% outputs at the current time instant given the data
% at the current time:
%
%   x_hat(k|k) : estimate of states at time k
%   y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states at the next time
% instant given the data at the current time are also
% calculated:
%
%   x_hat(k+1|k) : estimate of states at time k + 1
%
% The system model is defined as:
%
%   x(k+1) = A(k) x(k) + B(k) u(k) + w(k)
%     y(k) = C(k) x(k) + v(k)
%
% Note: there is no direct transmission (D = 0).
%
% The observer object can be used recursively in an 
% iteration loop or in a Simulink S-function block (see 
% MKFObserver_sfunc.m)
%
% Arguments:
%   models : (1, nj) cell array of structs
%       Each struct contains the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and Ts for the sample period.
%   P0 : (n, n) double
%       Initial covariance matrix of the state estimates
%       (same for each filter).
%   T : Transition probabity matrix of the Markov switching
%       process.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
%   p_seq_g_Yk_init : (nj, 1) (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%

classdef MKFObserverGPB2 < MKFObserver
    properties
        merged struct
    end
    methods
        function obj = MKFObserverGPB2(models,P0,T,label,x0, ...
                p_seq_g_Yk_init)
            arguments
                models (1, :) cell
                P0 double
                T double
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
            end

            % Get number of system models and check their dimensions
            nj = check_models(models);

            % Initial modes of system at time k = 0
            [~, r0] = mode_transitions_all(nj);
            % r0 is not used after initialization. Only needed here for
            % passing to MKFObserver to set the number of hypotheses

            if isempty(p_seq_g_Yk_init)
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(nj, 1) ./ nj;
            end

            % Split prior mode probabilities into nj^2 probabilities
            p_seq_g_Yk_init2 = reshape( ...
                repmat(p_seq_g_Yk_init ./ nj, nj, 1), ...
                [], 1 ...
            );

            % Create super-class observer instance
            obj = obj@MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init2,false);

            % Initialize merged estimates struct
            merged.p_seq_g_Yk_init = p_seq_g_Yk_init;

            % Add additional variables used by this observer
            obj.merged = merged;
            obj.type = "MKF_GPB2";

            % Initialize all variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Call reset method of super class
            reset@MKFObserver(obj);

            % Vectors of mode transitions to be modelled
            [obj.rkm1, obj.rk] = mode_transitions_all(obj.nj);

            % Create struct to store merged estimates and
            % initialize with initial state values
            obj.merged.Xk_est = nan(obj.n, 1, obj.nj);
            obj.merged.Pk = nan(obj.n, obj.n, obj.nj);
            obj.merged.Yk_est = nan(obj.ny, 1, obj.nj);
            obj.merged.p_seq_g_Yk = nan(obj.nj, 1);

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
        % updates the estimates of the switching Kalman filter
        % and calculates the predictions of the states and output
        % at the next sample time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %

            % Check size of arguments passed
            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

% External function based on code from Zhao et al. - NOT working yet
            % Update state and output estimates based on current
            % measurement and prior predictions
%             [obj.xk_est, obj.yk_est, obj.Pk, obj.merged.Xk_est, ...
%                 obj.merged.Yk_est, obj.merged.Pk, obj.merged.p_seq_g_Yk] = ...
%                 GPB2_update( ...
%                     obj.models, ...
%                     obj.T, ...
%                     obj.filters.Xkp1_est, ...
%                     obj.filters.Pkp1, ...
%                     yk, ...
%                     obj.p_seq_g_Yk ...
%                 );
%             assert(~any(isnan(obj.xk_est)))
%             assert(~any(isnan(obj.merged.Xk_est)))
%             assert(~any(isnan(obj.p_seq_g_Yk)))

            % Kalman filter update step
            % Inputs:
            %   - obj.filters.Xkp1_est
            %   - obj.filters.Pkp1
            %   - yk
            % Calculates:
            %   - obj.filters.Xk_est
            %   - obj.filters.Pk
            %   - obj.filters.Yk_est
            obj.KF_update(yk, obj.rk)

            % Calculate hypothesis probabilities (priors)
            obj.MKF_prob_prior(obj.rk, obj.rkm1)
            % Alternatively, could calculate obj.p_rk_g_rkm1 at
            % initialization and then just do this each time:
            %obj.p_seq_g_Ykm1 = obj.p_rk_g_rkm1 .* obj.p_seq_g_Yk;

            % Bayesian updating of hypothesis probabilities
            % Calculates:
            %   - obj.p_seq_g_Yk
            obj.MKF_prob_update(yk)
            assert(sum(obj.p_seq_g_Yk) - 1 < 1e-15)

            % Merge the nj^2 estimates into nj modes
            [obj.merged.Xk_est, obj.merged.Pk, obj.merged.Yk_est, ...
                obj.merged.p_seq_g_Yk] = ...
                merge_estimates( ...
                    obj.filters.Xk_est, ...
                    obj.filters.Pk, ...
                    obj.filters.Yk_est, ...
                    obj.p_seq_g_Yk, ...
                    obj.rk, ...
                    obj.rk ...
                );
            if any(isnan(obj.merged.Xk_est))
                disp("Stop")
            end
            %TODO: Consider removing these checks after testing
            assert(~any(isnan(obj.merged.Xk_est), 'all'))
            assert(~any(isnan(obj.merged.Pk), 'all'))
            assert(~any(isnan(obj.merged.Yk_est), 'all'))
            assert(~any(isnan(obj.merged.p_seq_g_Yk), 'all'))
            assert(sum(obj.merged.p_seq_g_Yk) - 1 < 1e-15)

            % Merge the estimates and error covariances again into 
            % a single set of estimates and error covariances
            % TODO: Or why not merge the original nj^2 estimates?
            [obj.xk_est, obj.Pk, obj.yk_est, p_check] = ...
                merge_estimates( ...
                    obj.merged.Xk_est, ...
                    obj.merged.Pk, ...
                    obj.merged.Yk_est, ...
                    obj.merged.p_seq_g_Yk ...
                );
            %TODO: Consider removing these checks after testing
            assert(~any(isnan(obj.xk_est), 'all'))
            assert(~any(isnan(obj.Pk), 'all'))
            assert(~any(isnan(obj.yk_est), 'all'))
            assert(sum(p_check) == 1)

            % GBP2 clones the nj merged estimates into nj^2 and 
            % uses these for making the nj^2 predictions at the
            % next time instant:
            %   xi_est(k+1|k) = Ai(k) * x_est(k|k-1) + Bi(k) * u(k);
            %   Pi(k+1|k) = Ai(k) * P(k|k-1) * Ai(k)' + Qi(k);
            obj.filters.Xk_est = obj.merged.Xk_est(obj.rkm1);
            obj.filters.Pk = obj.merged.Pk(obj.rkm1);

            % Kalman filter prediction step - computes estimates
            % of states and error covariances at the next time
            % instant
            % Inputs:
            %   - obj.filters.Xk_est
            %   - obj.filters.Pk
            %   - uk
            %   - rk
            % Calculates:
            %   - obj.filters.Xkp1_est
            %   - obj.filters.Pkp1
            obj.KF_predict(uk, obj.rk)

            % TODO: Do we still need merged xkp1_est and
            % Pkp1 estimates?
            weights = reshape(obj.p_seq_g_Yk, 1, 1, []);
            obj.xkp1_est = sum(weights .* obj.filters.Xkp1_est, 3);
            assert(~any(isnan(obj.xkp1_est)))
            Xkp1_devs = obj.xkp1_est - obj.filters.Xkp1_est;
            obj.Pkp1 = sum(weights .* (obj.filters.Pkp1 + ...
                pagemtimes(Xkp1_devs, pagetranspose(Xkp1_devs))), 3);

            %TODO: Should the probabilities be cloned from merged
            % probabilities? or leave the unmerged p_seq_g_Yk values 
            % from current period. Yes. See eq.n 37 in Zhou et al.
            obj.p_seq_g_Yk = obj.merged.p_seq_g_Yk(obj.rkm1);
            %TODO: Should these be normalized?  Don't think so
            % because in next time instant, obj.p_seq_g_Yk will
            % be scaled when it is updated.  Not necessary, makes
            % no difference.
            %obj.p_seq_g_Yk = obj.p_seq_g_Yk ./ sum(obj.p_seq_g_Yk);

        end
    end
end