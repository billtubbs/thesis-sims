% Multi-model Kalman Filter class definition
%
% obs = MKFObserverSF_DI(models,P0,seq,T,d,label,x0,r0, ...
%     p_seq_g_Yk_init,reset)
%
% Class for simulating a sub-optimal multi-model observer with 
% sequence fusion, for state estimation of a Markov jump linear 
% system. This version includes the detection interval 
% procedure described by Robertson et al. (1998) which reduces 
% the frequency of the branching and pruning procedure. 
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next
% time instant are also calculated:
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
%   seq : model indicator sequences for each filter (in rows).
%   T : Transition probabity matrix for the Markov switching
%       process. T(i,j) represents the probability of the
%       system switching from model i to j.
%   d : integer double
%      Detection interval length in number of sample periods.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double
%       Initial state estimates (optional, default zeros).
%   r0 : (1, 1) or (nh, 1) integer (optional)
%       Integer scalar or vector with values in the range 
%       {1, ..., nj} which indicate the prior system modes at time 
%       k = -1. If not provided, the default initialization based 
%       on the mode sequence that is generated will be used.
%   p_seq_g_Yk_init : (nh, 1) double (optional, default uniform)
%       Initial prior probabilities of each hypothesis at
%       time k-1. If not specified, default is equal, i.e.
%       uniform, probability assigned to each hypothesis.
%

classdef MKFObserverSF_DI < MKFObserverSF
    properties
        f double {mustBeInteger, mustBeNonnegative}
        d double {mustBeInteger, mustBeNonnegative}
        id int16
        id_next int16
    end
    methods
        function obj = MKFObserverSF_DI(models,P0,seq,T,d,label,x0, ...
            r0,p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                seq (:, 1) cell
                T double
                d (1, 1) double {mustBeInteger, ...
                    mustBeGreaterThanOrEqual(d,1)}
                label (1, 1) string = ""
                x0 = []
                r0 (:, 1) double {mustBeInteger, ...
                    mustBeGreaterThanOrEqual(r0, 1)} = []
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Create MKF super-class observer instance
            obj = obj@MKFObserverSF(models,P0,seq,T,label,x0,r0, ...
                p_seq_g_Yk_init,false)

            % Add additional variables used by SP observer
            obj.d = d;
            obj.f = obj.nf .* d;
            obj.type = "MKF_SF_DI";

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

            % Call reset method of super class object
            reset@MKFObserverSF(obj);

            % Reset sequence index
            obj.id = int16(0);
            obj.id_next = int16(1);

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
        % updates the estimates of the multi-model Kalman filter
        % and calculates the predictions of the states and output
        % at the next sample time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %

            % Increment counter and sequence index:
            % obj.id is the counter for the detection interval 
            % (1 <= id <= d). Switching is assumed to only 
            % occur at the end of the detection interval.
            % When obj.id exceeds d (the spacing parameter), it 
            % is reset to 1.
            obj.id = obj.id_next;
            obj.id_next = mod(obj.id, obj.d) + 1;
            if obj.id == 1
                % Increment sequence index (at end of sequence it 
                % wraps to beginnning)
                obj.i = obj.i_next;
                obj.i_next = mod(obj.i, obj.nf) + 1;
            end

            % If at end of a detection interval, carry out
            % hypothesis merging and branching procedure and 
            % update conditional probability estimates.
            if obj.id == obj.d

                % Run MKF super-class updates and probability calcs
                %update@MKFObserverSF(obj, yk, uk);

%                 update@MKFObserverBM(obj, yk, uk, ...
%                     obj.idx_modes{obj.i}, ...
%                     obj.idx_merge{obj.i}, ...
%                     obj.idx_branch{obj.i_next});

                % Instead of calling above function, implement it here
                % Copy these as if the function had been called
                idx_modes = obj.idx_modes{obj.i};
                idx_merge = obj.idx_merge{obj.i};
                idx_branch = obj.idx_branch{obj.i_next};

                % Set r(k-1) to previous values of r(k)
                rkm1 = obj.rk;
    
                % Current modes from sequence index
                rk = idx_modes;
    
                % Re-initialize arrays
                obj.filters.Yk_est = nan(obj.ny, 1, obj.nh);
                obj.filters.Kf = nan(obj.n, obj.ny, obj.nh);
                obj.filters.Sk = nan(obj.ny, obj.ny, obj.nh);
                obj.p_yk_g_seq_Ykm1 = nan(obj.nh, 1);

                % Kalman filter update step
                obj.KF_update(yk, rk)

                % TODO: WHOLE REASON FOR COPYING THIS CODE HERE
                % WAS SIMPLY TO SKIP THIS STEP. THERE MUST BE A
                % AN EASIER WAY!
                % Calculate hypothesis probabilities (priors)
                %obj.MKF_prob_prior(rk, rkm1)
                % Instead set prior to posterior probabilities
                % calculated at last iteration (equivalent to
                % saying transition probabilities are all 1)
                obj.p_seq_g_Ykm1 = obj.p_seq_g_Yk;

                % Bayesian updating of hypothesis probabilities
                obj.MKF_prob_update(yk, rk)

                % TODO: Not sure how to do this as error obs> 5e-15 can happen
                assert(sum(obj.p_seq_g_Yk) - 1 < 1e-12)

                % Partialy merge the nh estimates and error covariances
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

                % Fully merge the estimates and error covariances
                [obj.xk_est, obj.Pk, obj.yk_est, p_check] = ...
                    merge_estimates( ...
                        obj.merged.Xk_est, ...
                        obj.merged.Pk, ...
                        obj.merged.Yk_est, ...
                        obj.merged.p_seq_g_Yk ...
                    );

                % Branch (i.e. clone) the nm merged estimates
                obj.filters.Xk_est = obj.merged.Xk_est(:,:,idx_branch);
                obj.filters.Pk = obj.merged.Pk(:,:,idx_branch);
                obj.p_seq_g_Yk = obj.merged.p_seq_g_Yk(idx_branch);
                rk = obj.merged.rk(idx_branch);
                %TODO: Should the prob's be normalized here?
                obj.p_seq_g_Yk = obj.p_seq_g_Yk ./ sum(obj.p_seq_g_Yk);

                % Update number of branched hypothesis (may be time-varying)
                if obj.nh ~= length(rk)
                    obj.nh = length(rk);
                end

                % Kalman filter prediction step
                obj.KF_predict(uk, rk)

                % Store mode transitions 
                obj.rkm1 = rkm1;
                obj.rk = rk;

            else

                rk = obj.idx_modes{obj.i};

                % At every other sample time within the detection 
                % interval, simply update Kalman filters and merged 
                % estimates assuming the same (i.e. repeating) mode 
                % transitions
                if obj.id == 1
                    % At start of detection interval, set r(k-1) to
                    % the modes from the previous detection interval.
                    rkm1 = obj.rk;
                    % Compute priors - to be used throughout detection
                    % interval
                    obj.MKF_prob_prior(rk, rkm1)
                else
                    % Otherwise, reset to same modes
                    rkm1 = obj.rkm1;
                    % and reset prior to posterior probabilities
                    % calculated at last iteration (equivalent to
                    % saying transition probabilities are all 1)
                    obj.p_seq_g_Ykm1 = obj.p_seq_g_Yk;
                end

                %update@MKFObserver(obj, yk, uk, rk, rkm1)

                % Carry out all steps except prior update
                obj.KF_update(yk, rk)

                % Don't recompute prior since there is no
                % mode transition
                %obj.MKF_prob_prior(rk, rkm1)

                obj.MKF_prob_update(yk, rk)
                obj.MKF_estimates()
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
                obj.rkm1 = rkm1;
                obj.rk = rk;

                % This is only needed for producing outputs for 
                % plotting
                idx_merge = obj.idx_merge{obj.i};
                % Merge the nh estimates into nm - 
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

            end

        end
    end
end
