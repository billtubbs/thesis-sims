% Multi-model Kalman Filter class definition
%
% obs = MKFObserverBM(models,P0,T,r0,label,x0, ...
%     p_seq_g_Yk_init,reset)
%
% Class for simulating a multi-model Kalman filter for 
% state estimation of a Markov jump linear system. 
% 
% This version is based on MKFObserver but provides 
% branching and merging methods used by sub-optimal 
% algorithms.
%
% This is the filtering form of the observer, which 
% produces posterior estimates of the states and outputs 
% at the current time instant given the data at the 
% current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next
% time instant given the data at the current time are
% also calculated:
%
%  x_hat(k+1|k) : estimate of states at time k + 1
%  y_hat(k+1|k) : estimate of outputs at time k + 1
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
%   T : Transition probabity matrix for the Markov switching
%       process. T(i,j) represents the probability of the
%       system switching from model i to j.
%   r0 : (nh, 1) integer (optional, default zeros)
%       Initial prior model indicator value at time k-1.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double
%       Initial state estimates (optional, default zeros).
%   p_seq_g_Yk_init : (nh, 1) double (optional, default uniform)
%       Initial prior probabilities of each hypothesis at
%       time k-1. If not specified, default is equal, i.e.
%       uniform, probability assigned to each hypothesis.
%   reset : logical (default, true)
%       If true, the object's reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       reseting).
%

classdef MKFObserverBM < MKFObserver
    properties
        merged struct
        nm double
    end
    methods
        function obj = MKFObserverBM(models,P0,T,r0,label,x0, ...
            p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                T double
                r0 (:, 1) int16 {mustBeGreaterThan(r0, 0)}
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Create MKF super-class observer instance
            obj = obj@MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,false)

            % Add additional variables used by this observer
            obj.merged = struct();
            obj.type = "MKF_BM";

            if reset
                % Initialize variables
                obj.reset()
            end

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Call reset method of super class
            reset@MKFObserver(obj);

            % Create struct to store merged estimates and
            % covariances. Note: the size of these may change
            % each time step if nm is not the same.
            obj.merged.Xk_est = [];  % nan(obj.n, 1, obj.nm);
            obj.merged.Pk = [];  % nan(obj.n, obj.n, obj.nm);
            obj.merged.Yk_est = [];  % nan(obj.ny, 1, obj.nm);
            obj.merged.p_seq_g_Yk = [];  % nan(obj.nm, 1);
            obj.merged.rk = [];  % nan(obj.nm, 1);

        end
        function update(obj, yk, uk, idx_modes, idx_merge, idx_branch)
        % obj.update(yk, uk, idx_merge, idx_branch, idx_modes)
        % carries out the hypothesis merging, Kalman filter 
        % update step, hypothesis branching and Kalman filter
        % prediction steps and calculates the mixed estimates 
        % of the states and outputs at the current sample time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %   idx_merge : vector (nh, 1) of hypotheses indeces to
        %       use for hypothesis merging.
        %   idx_branch : vector (nh, 1) of hypotheses indeces to
        %       use for hypothesis branching.
        %   idx_modes : vector, size (nj, 1)
        %       System modes at current time k.
        %
        % Note the order of merging, branching, and mode transitions 
        % used here:
        %  - idx_modes are the modes to use for the KF prediction
        %    and update steps in the current time step to calculate 
        %    x_est_f(k|k), y_est(k|k), and P(k|k), and also to 
        %    calculate the transition probabilities.
        %  - idx_merge is then used to merge the updated estimates and 
        %    covariances, using the conditoinal probabilities, 
        %    Pr(Gamma_f(k)|Y(k)).
        %  - Finally, idx_branch is used to branch the merged 
        %    hypotheses before the KF prediction step. However, it must 
        %    match the modes and merging indeces that will occur in the 
        %    next time step. Therefore, idx_branch(k+1) should be 
        %    provided here, not idx_branch(k).
        %

            % Check size of arguments passed
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")
            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(idx_modes), [obj.nh 1]), "ValueError: size(idx_modes)")
            assert(isequal(size(idx_merge), [obj.nh 1]), "ValueError: size(idx_merge)")

            % Set r(k-1) to previous values of r(k)
            rkm1 = obj.rk;

            % Current modes from sequence index
            rk = idx_modes;

            % Since nh is possibly time-varying, the array sizes
            % may have changed so re-initialize before calling
            % KF_update and MKF_prob_update
            obj.filters.Yk_est = nan(obj.ny, 1, obj.nh);
            obj.filters.Kf = nan(obj.n, obj.ny, obj.nh);
            obj.filters.Sk = nan(obj.ny, obj.ny, obj.nh);
            obj.p_yk_g_seq_Ykm1 = nan(obj.nh, 1);
            % (Note that obj.p_seq_g_Ykm1 and obj.p_seq_g_Yk 
            %  will change size automatically when
            %  MKF_prob_update so this is not necessary)

            % Kalman filter update step
            % Inputs:
            %   - obj.filters.Xkp1_est
            %   - obj.filters.Pkp1
            %   - yk
            % Calculates:
            %   - obj.filters.Xk_est
            %   - obj.filters.Pk
            %   - obj.filters.Yk_est
            obj.KF_update(yk, rk)

            % TODO: Could this method simply call the parent class method?
            % Run MKF super-class updates and probability calcs
            % update@MKFObserver(obj, yk, uk, rk);

            % Calculate hypothesis probabilities (priors)
            obj.MKF_prob_prior(rk, rkm1)

            % Bayesian updating of hypothesis probabilities
            % Calculates:
            %   - obj.p_seq_g_Yk
            obj.MKF_prob_update(yk, rk)
            % TODO: Not sure how to do this as error obs> 5e-15 can happen
            assert(sum(obj.p_seq_g_Yk) - 1 < 1e-12)

            % Merge the nh estimates into nm
            [obj.merged.Xk_est, obj.merged.Pk, obj.merged.Yk_est, ...
                obj.merged.p_seq_g_Yk, obj.merged.rk] = ...
                merge_estimates( ...
                    obj.filters.Xk_est, ...
                    obj.filters.Pk, ...
                    obj.filters.Yk_est, ...
                    obj.p_seq_g_Yk, ...
                    rk, ...
                    idx_merge ...
                );
            %TODO: Consider removing these checks after testing
            %assert(~any(isnan(obj.merged.Xk_est), 'all'))
            %assert(~any(isnan(obj.merged.Pk), 'all'))
            %assert(~any(isnan(obj.merged.Yk_est), 'all'))
            %assert(~any(isnan(obj.merged.p_seq_g_Yk), 'all'))
            %assert(sum(obj.merged.p_seq_g_Yk) - 1 < 1e-6)
            %assert(isequal(obj.merged.rk, ...
            %    cellfun(@(s) s(obj.i), obj.seq)), ...
            %    "ValueError: idx_merge not consistent with seq")

            % Merge the estimates and error covariances again into 
            % a single set of estimates and error covariances
            [obj.xk_est, obj.Pk, obj.yk_est, p_check] = ...
                merge_estimates( ...
                    obj.merged.Xk_est, ...
                    obj.merged.Pk, ...
                    obj.merged.Yk_est, ...
                    obj.merged.p_seq_g_Yk ...
                );
            % TODO: Or why not merge the original nj^2 estimates?
            %[obj.xk_est, obj.Pk, obj.yk_est, ~] = ...
            %    merge_estimates( ...
            %        obj.filters.Xk_est, ...
            %        obj.filters.Pk, ...
            %        obj.filters.Yk_est, ...
            %        obj.p_seq_g_Yk ...
            %    );

            %TODO: Consider removing these checks after testing
            %assert(isequal(size(obj.xk_est), [obj.n 1]))
            %assert(isequal(size(obj.Pk), [obj.n obj.n]))
            %assert(isequal(size(obj.yk_est), [obj.ny 1]))
            %assert(~any(isnan(obj.xk_est), 'all'))
            %assert(~any(isnan(obj.Pk), 'all'))
            %assert(~any(isnan(obj.yk_est), 'all'))
            %assert(p_check == 1)

            % Branch (i.e. clone) the nm merged estimates and mode
            % index into nh estimates and modes to be used for making 
            % the nh predictions at the next time instant:
            %   xi_est(k+1|k) = Ai(k) * x_est(k|k-1) + Bi(k) * u(k);
            %   Pi(k+1|k) = Ai(k) * P(k|k-1) * Ai(k)' + Qi(k);
            obj.filters.Xk_est = obj.merged.Xk_est(:,:,idx_branch);
            obj.filters.Pk = obj.merged.Pk(:,:,idx_branch);
            obj.p_seq_g_Yk = obj.merged.p_seq_g_Yk(idx_branch);
            rk = obj.merged.rk(idx_branch);
            %TODO: Should the prob's be normalized here?
            obj.p_seq_g_Yk = obj.p_seq_g_Yk ./ sum(obj.p_seq_g_Yk);

            % Update number of branched hypothesis (may be time-varying)
            if obj.nh ~= length(rk)
                obj.nh = length(rk);
                % TODO: Should we delete these estimates as they will 
                % no longer be valid after branching and need to be 
                % re-calculated in the next time step?
                % No: leave them for testing/debugging and reset their
                % size in MKFObserver.update at next time step.
                %obj.filters.Yk_est = nan(obj.ny, 1, obj.nh);
                %obj.p_yk_g_seq_Ykm1 = nan(obj.nh, 1);
            end

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
            obj.KF_predict(uk, rk)

            % TODO: Do we still need merged xkp1_est and
            % Pkp1 estimates? Delete. These cannot be merged yet
            % because the mode transitions are not known until
            % the next time step.
            %weights = reshape(obj.p_seq_g_Yk, 1, 1, []);
            %obj.xkp1_est = sum(weights .* obj.filters.Xkp1_est, 3);
            %assert(~any(isnan(obj.xkp1_est)))
            %Xkp1_devs = obj.xkp1_est - obj.filters.Xkp1_est;
            %obj.Pkp1 = sum(weights .* (obj.filters.Pkp1 + ...
            %    pagemtimes(Xkp1_devs, pagetranspose(Xkp1_devs))), 3);

            % Store mode transitions (obj.rk is used in the next 
            % iteration). However, note that rkm1 here are the modes  
            % before the merging and branching steps that occurred at 
            % the current time and therefore rk does not correspond
            % to the transitions from rkm1. This is okay since rkm1 
            % is not needed and is only stored for debugging purposes.
            % (actually, MKFObserverSF_RODD uses it to reset the
            %  transitions when repeating them during the detection
            %  interval).
            obj.rkm1 = rkm1;
            obj.rk = rk;

        end

    end
end
