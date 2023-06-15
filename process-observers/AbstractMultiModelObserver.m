% Abstract class definition for multiple-model filters
%
% obs = AbstractMultiModelObserver(models,type,T,r0,label, ...
%    p_seq_g_Yk_init,reset)
%
% Abstract Class for inheriting when defining multiple-model 
% observers for state estimation of Markov jump linear systems. 
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
% Arguments:
%   models : (1, nj) cell array of structs
%       Each struct contains the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and Ts for the sample period.
%   type : string
%       Code used to categorize the specific type of 
%       filter.
%   T : Transition probabity matrix of the Markov switching
%       process.
%   r0 : (nh, 1) integer
%       Integer vector with values in the range {1, ..., nj} 
%       which indicates the prior system mode at time k = -1.
%   label : string
%       Arbitrary name to identify the observer.
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

classdef (Abstract) AbstractMultiModelObserver < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        nh (1, 1) double {mustBeInteger, mustBeNonnegative}
        models (1, :) cell
        Ts (1, 1) double {mustBeNonnegative}
        T double
        r0 (:, 1) int16 {mustBeGreaterThanOrEqual(r0, 1)}
        label (1, 1) string
        p_seq_g_Yk_init double
        p_seq_g_Ykm1 double
        p_seq_g_Yk double
        p_yk_g_seq_Ykm1 double
        p_rk_g_Ykm1 double
        p_rk_g_rkm1 double
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        Pkp1 double
        rk (:, 1) int16
        rkm1 (:, 1) int16
        type (1, 1) string
    end
    methods
        function obj = AbstractMultiModelObserver(models,type,T,r0,label, ...
                p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                type (1, 1) string
                T double
                r0 (:, 1) double {mustBeInteger, ...
                    mustBeGreaterThanOrEqual(r0, 1)}
                label (1, 1) string = ""
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Number of models
            nj = length(models);

            % Number of hypotheses to be modelled
            nh = size(r0, 1);

            if label == ""
                label = type;
            end
            if isempty(p_seq_g_Yk_init)
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(nh, 1) ./ nh;
            end

            % Check transition probability matrix
            assert(isequal(size(T), [nj nj]), "ValueError: size(T)")
            assert(all(abs(sum(obj.T, 2) - 1) < 1e-15), "ValueError: T")

            % Store parameters
            obj.T = T;
            obj.nj = nj;
            obj.models = models;
            obj.nh = nh;
            obj.r0 = int16(r0);
            obj.p_seq_g_Yk_init = p_seq_g_Yk_init;
            obj.type = type;
            obj.label = label;

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

            % Initialize state and output estimates
            obj.xkp1_est = obj.x0;
            obj.Pkp1 = obj.P0;
            obj.rk = obj.r0;
            obj.rkm1 = zeros(size(obj.r0));  % nan not possible for int16

            % Number of hypotheses - may be time-varying
            obj.nh = size(obj.rk, 1);

            % Initial values of prior conditional probabilities at k = -1 
            obj.p_seq_g_Yk = obj.p_seq_g_Yk_init;

            % Empty vectors to store values for filter calculations
            % p(y(k)|R(k),Y(k-1))
            obj.p_yk_g_seq_Ykm1 = nan(obj.nh, 1);
            % Pr(r(k)|Y(k-1))
            obj.p_rk_g_Ykm1 = nan(obj.nh, 1);
            % Pr(R(k))
            obj.p_rk_g_rkm1 = nan(obj.nh, 1);
            % Pr(R(k)|Y(k-1))
            obj.p_seq_g_Ykm1 = nan(obj.nh, 1);

            % At initialization at time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);

        end
        function MKF_prob_prior(obj, rk, rkm1)
        % Calculate prior probabilities of current mode 
        % transitions

            % Compute Pr(r(k)|r(k-1)) based on Markov transition
            % probability matrix
            obj.p_rk_g_rkm1 = prob_rk(rk, obj.T(rkm1, :)');

            % Compute Pr(r(k)|Y(k-1)) in current timestep from
            % previous estimate (Pr(r(k-1)|Y(k-1))) and transition
            % probabilities
            obj.p_seq_g_Ykm1 = obj.p_rk_g_rkm1 .* obj.p_seq_g_Yk;

        end
        function MKF_prob_update(obj, yk, rk)
        % Bayesian updates to conditional probability
        % estimates of each hypothesis

            if nargin < 3
                rk = obj.rk;
            end

            % Estimate conditional probabilities of predictions
            % for each hypothesis
            for f = 1:obj.nh

                % Select system model for this mode
                m = obj.models{rk(f)};

                % Output prediction y_f_est(k|k-1)
                yk_pred = m.C * obj.filters.Xkp1_est(:,:,f);

                % Covariance of the output prediction error
                cov_yk = m.C * obj.filters.Pkp1(:,:,f) * m.C' + m.R;

                % Ensure covariance matrix is symmetric
                if ~isscalar(cov_yk)
                    cov_yk = (cov_yk + cov_yk.')/2;
                end

%                 % Check if it is symmteric, positive definite
%                 tf = issymmetric(cov_yk);
%                 if ~tf
%                     error("ValueError: cov_yk not symmetric")
%                 end
%                 d = eig(cov_yk);
%                 isposdef = all(d > 0);
%                 if ~isposdef
%                     error("ValueError: cov_yk not positive definite")
%                 end

                % Probability density of output prediction (assuming a
                % multivariate normal distribution) based on prior
                % estimates computed in previous timestep
                obj.p_yk_g_seq_Ykm1(f) = mvnpdf(yk, yk_pred, cov_yk);

            end
            assert(~any(isnan(obj.p_yk_g_seq_Ykm1)))
            assert(~all(obj.p_yk_g_seq_Ykm1 == 0))

            % Bayesian update to calculate posterior, Pr(r(k)|Y(k))
            cond_pds = obj.p_yk_g_seq_Ykm1 .* obj.p_seq_g_Ykm1;
            obj.p_seq_g_Yk = cond_pds ./ sum(cond_pds);
            % Note: above calculation normalizes p_seq_g_Yk so that
            % assert(abs(sum(obj.p_seq_g_Yk) - 1) < 1e-15) % is always true
            %fprintf("[%d %d]: MKF_prob_update [%d %d %d %d %d %d] -> [%d %d %d %d %d %d]\n",  obj.i,  obj.id, rkm1, rk)

        end
        function MKF_estimates(obj)
        % Merge the Kalman filter estimates and error covariances
        % using a weighted-average based on the hypothesis
        % conditional probabilities.
            [obj.xk_est, obj.Pk, obj.yk_est, p_check] = ...
                merge_estimates( ...
                    obj.filters.Xk_est, ...
                    obj.filters.Pk, ...
                    obj.filters.Yk_est, ...
                    obj.p_seq_g_Yk ...
                );
            assert(p_check == 1)
        end
    end
end
