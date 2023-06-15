% Multi-model Kalman Filter class definition
%
% obs = MKFObserver(models,P0,T,r0,label,x0,p_seq_g_Yk_init,reset)
%
% Class for simulating a multi-model Kalman filter for state
% estimation of a Markov jump linear system. 
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
%   T : Transition probabity matrix of the Markov switching
%       process.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
%   r0 : (nh, 1) integer
%       Integer vector with values in the range {1, ..., nj} 
%       which indicates the prior system mode at time k = -1.
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%   reset : logical (default, true)
%       If true, the object's reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       reseting).
%

classdef MKFObserver < AbstractMultiModelObserver
    properties (SetAccess = immutable)
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        P0 double
        x0 (:, 1) double
        filters  struct
    end
    methods
        function obj = MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                T double
                r0 (:, 1) double {mustBeInteger, ...
                    mustBeGreaterThanOrEqual(r0, 1)}
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Get number of system models and check their dimensions
            [nj, n, nu, ny, Ts] = check_models(models);

            % Call super-class constructor
            obj = obj@AbstractMultiModelObserver(models,"MKF",T,r0, ...
                label,p_seq_g_Yk_init,false);

            % Check dimensions of other parameters
            assert(isequal(size(P0), [n n]), "ValueError: size(P0)")
            for j = 1:nj
                assert(isequal(size(models{j}.Q), [n n]))
                assert(isequal(size(models{j}.R), [ny ny]))
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]), "ValueError: x0")
            end

            % Create struct to store Kalman filter variables
            obj.filters = struct();
            obj.filters.Xkp1_est = nan(n, 1, obj.nh);
            obj.filters.Pkp1 = nan(n, n, obj.nh);
            obj.filters.Xk_est = nan(n, 1, obj.nh);
            obj.filters.Pk = nan(n, n, obj.nh);
            obj.filters.Yk_est = nan(ny, 1, obj.nh);
            obj.filters.Kf = nan(n, ny, obj.nh);
            obj.filters.Sk = nan(ny, ny, obj.nh);

            % Store parameters
            obj.Ts = Ts;
            obj.nu = nu;
            obj.n = n;
            obj.ny = ny;
            obj.P0 = P0;
            obj.x0 = x0;

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

            % Call super-class reset method
            reset@AbstractMultiModelObserver(obj)

            % Reset Kalman filter variables
            % Note: At initialization at time k = 0, Xkp1_est and
            % Pkp1 represent prior estimates of the states, and
            % covariance, i.e. x_est_f(k|k-1) and P_f(k|k-1).
            obj.filters.Xkp1_est = repmat(obj.x0, 1, 1, obj.nh);
            obj.filters.Pkp1 = repmat(obj.P0, 1, 1, obj.nh);
            obj.filters.Xk_est(:,:,:) = nan;
            obj.filters.Pk(:,:,:) = nan;
            obj.filters.Yk_est(:,:,:) = nan;
            obj.filters.Kf(:,:,:) = nan;
            obj.filters.Sk(:,:,:) = nan;

        end
        function KF_update(obj, yk, rk)
        % Update Kalman filter estimates using current 
        % measurement, y(k).  Also stores the current values
        % of the correction gains and output error covariances
        % in obj.filters.Kf and obj.filters.Sk for debugging
        % and analysis.
            for f = 1:obj.nh
                m = obj.models{rk(f)};
                [obj.filters.Xk_est(:,:,f), obj.filters.Pk(:,:,f), ...
                    obj.filters.Yk_est(:,:,f), obj.filters.Kf(:,:,f), ...
                    obj.filters.Sk(:,:,f)] = ...
                        kalman_update_f( ...
                            m.C, m.R, ...
                            obj.filters.Xkp1_est(:,:,f), ...
                            obj.filters.Pkp1(:,:,f), ...
                            yk ...
                        );
            end
            % TODO: For testing only - remove later
            %assert(~any(isnan(obj.filters.Xk_est), 'all'))
            %assert(~any(isnan(obj.filters.Pk), 'all'))
        end
        function KF_predict(obj, uk, rk)
        % Calculate Kalman filter predictions of
        % states at next time instant using current
        % input u(k).
            for f = 1:obj.nh
                m = obj.models{rk(f)};
                [obj.filters.Xkp1_est(:,:,f), obj.filters.Pkp1(:,:,f)] = ...
                    kalman_predict_f( ...
                        m.A, m.B, m.Q, ...
                        obj.filters.Xk_est(:,:,f), ...
                        obj.filters.Pk(:,:,f), ...
                        uk ...
                    );
            end
            % TODO: For testing only - remove later
            %assert(~any(isnan(obj.filters.Xkp1_est(:,:,f)), 'all'))
            %assert(~any(isnan(obj.filters.Pkp1(:,:,f)), 'all'))
        end
        function update(obj, yk, uk, rk, rkm1)
        % obj.update(yk, uk, rk, rkm1, step)
        % updates the estimates of the multi-model Kalman filter
        % and calculates the predictions of the states and outputs
        % at the next sample time.
        %
        % Arguments:
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %   rk : vector, size (nj, 1)
        %       System modes at current time k.
        %   rkm1 : vector, size (nj, 1) (optional)
        %       System modes at time k - 1. If not specified,
        %       rkm1 is set to the values stored in obj.rk
        %       from the last call to this function.
        %

            % Check size of arguments passed
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")
            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(rk), [obj.nh 1]), "ValueError: size(rk)")

            if nargin < 5
                % Default: set r(k-1) to previous values of r(k)
                rkm1 = obj.rk;
            end

            % Kalman filter update step
            obj.KF_update(yk, rk)

            % Calculate hypothesis probabilities (priors)
            obj.MKF_prob_prior(rk, rkm1)

            % Bayesian updating of hypothesis probabilities
            obj.MKF_prob_update(yk, rk)
            % TODO: Remove this test
            assert(sum(obj.p_seq_g_Yk) - 1 < 5e-15)

            % Merge estimates and error covariances based on the 
            % hypothesis probabilities. 
            % TODO: Set this up so that sub-optimal algorithms
            % may override this method to also carry out hypothesis
            % merging, pruning and branching procedures, etc. 
            obj.MKF_estimates()

            % Kalman filter prediction step - computes estimates
            % of states and error covariances at the next time
            % instant
            obj.KF_predict(uk, rk)

            % TODO: Do we still need a merged xkp1 estimate?
            % I think not. It is meaningless in the cases where
            % merging and branching has occurred—unless the  
            % probabilities of branches are divided perhaps?
            [obj.xkp1_est, obj.Pkp1, ~, ~] = ...
                merge_estimates( ...
                    obj.filters.Xkp1_est, ...
                    obj.filters.Pkp1, ...
                    [], ...
                    obj.p_seq_g_Yk ...
                );

            %weights = reshape(obj.p_seq_g_Yk, 1, 1, []);
            %obj.xkp1_est = sum(weights .* obj.filters.Xkp1_est, 3);
            %assert(~any(isnan(obj.xkp1_est)))
            %Xkp1_devs = obj.xkp1_est - obj.filters.Xkp1_est;
            %obj.Pkp1 = sum(weights .* (obj.filters.Pkp1 + ...
            %    pagemtimes(Xkp1_devs, pagetranspose(Xkp1_devs))), 3);

            % Store branched mode transitions, r(k-1) -> r(k) 
            % that were used in the current time step (obj.rk
            % is used in the next timestep as r(k-1)).
            obj.rkm1 = rkm1;
            obj.rk = rk;

        end
    end
end
