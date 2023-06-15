% Multi-model Kalman Filter class definition
%
% obs = MKFObserverF(models,P0,T,r0,label,x0,p_seq_g_Yk_init)
% Class for simulating a multi-model Kalman filter for state
% estimation of a Markov jump linear system. This version
% differs from MKFObserver in that it uses instances of
% KalmanFilterF to do the hypothesis filter updates and
% predictions.
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
%   P0 : (n, n) double
%       Initial covariance matrix of the state estimates
%       (same for each filter).
%   T : Transition probabity matrix of the Markov switching
%       process.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
%   r0 : (nh, 1) integer (optional, default ones)
%       Integer in the range {1, ..., nj} which indicates
%       the prior system mode at time k = -1.
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%   reset : logical (default, true)
%       If true, the objects reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       reseting).
%

classdef MKFObserverF < MKFObserver
    methods
        function obj = MKFObserverF(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                T double
                r0 (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(r0, 1)}
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Get number of system models and check their dimensions
            [nj, n, nu, ny, Ts] = check_models(models);

            % Number of hypotheses to be modelled
            nh = size(r0, 1);

            % Check dimensions of other parameters
            assert(isequal(size(P0), [n n]), "ValueError: size(P0)")
            for j = 1:nj
                assert(isequal(size(models{j}.Q), [n n]))
                assert(isequal(size(models{j}.R), [ny ny]))
            end

            if isempty(p_seq_g_Yk_init)
                % Initial values of prior conditional probabilities at 
                % k = -1. In absence of prior knowledge, assume all 
                % equally likely
                p_seq_g_Yk_init = ones(nh, 1) ./ nh;
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]), "ValueError: x0")
            end

            % Create MKF super-class observer instance
            obj = obj@MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,false);

            % Create struct to store Kalman filter variables
            obj.filters = struct();
            obj.filters.Xkp1_est = nan(n, 1, obj.nh);
            obj.filters.Pkp1 = nan(n, n, obj.nh);
            obj.filters.Xk_est = nan(obj.n, 1, obj.nh);
            obj.filters.Pk = nan(obj.n, obj.n, obj.nh);
            obj.filters.Yk_est = nan(obj.ny, 1, obj.nh);

            % Modify selected parameters
            obj.type = "MKF_F";
            if label == ""
                label = obj.type;
            end
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

            % Initialize estimate covariance
            obj.Pkp1 = obj.P0;

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.Pkp1 = obj.P0;
            obj.rk = obj.r0;
            obj.rkm1 = [];

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

            % Reset Kalman filter variables
            obj.filters.Xkp1_est = repmat(obj.xkp1_est, 1, 1, obj.nh);
            obj.filters.Pkp1 = repmat(obj.Pkp1, 1, 1, obj.nh);
            obj.filters.Xk_est = nan(obj.n, 1, obj.nh);
            obj.filters.Pk = nan(obj.n, obj.n, obj.nh);
            obj.filters.Yk_est = nan(obj.ny, 1, obj.nh);

            % Add a cell array of switching Kalman Filters
            obj.filters.objects = cell(obj.nh, 1);
            fmt = strcat('%s%0', ...
                char(string(strlength(string(obj.nh)))), 'd');
            for i = 1:obj.nh
                label = sprintf(fmt,obj.label,i);
                % Index of system model
                obj.filters.objects{i} = ...
                    SKFObserver( ...
                        obj.models, ...
                        obj.P0, ...
                        label, ...
                        obj.x0, ...
                        obj.rk(i) ...
                    );
            end

            % Initialize state and output estimates
            % Note: At initialization at time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.Pkp1 = obj.P0;
            obj.rk = obj.r0;
            obj.rkm1 = [];

            % At initialization at time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);

        end
        function KF_update(obj, yk, rk)
        % Update Kalman filter estimates using current 
        % measurement, y(k).
            for f = 1:obj.nh
                kf = obj.filters.objects{f};

                % Update system mode, r(k)
                kf.rk = rk(f);

                % System model
                m = kf.models{rk(f)};

                % Update estimates based on current measurement
                [kf.xk_est, kf.Pk, kf.yk_est, kf.Sk, kf.Kf] = ...
                    kalman_update_f(m.C, ...
                    m.R, kf.xkp1_est, kf.Pkp1, yk);
                % TODO: For testing only - remove later
                assert(~any(isnan(kf.xk_est), 'all'))
                assert(~any(isnan(kf.Pk), 'all'))

                % Copy updated estimates to array
                obj.filters.Xk_est(:,:,f) = kf.xk_est;
                obj.filters.Pk(:,:,f) = kf.Pk;
                obj.filters.Yk_est(:,:,f) = kf.yk_est;

            end
        end
        function KF_predict(obj, uk, rk)
        % Calculate Kalman filter predictions of
        % states at next time instant using current
        % input u(k).
            for f = 1:obj.nh
                kf = obj.filters.objects{f};

                % Update system mode, r(k)
                kf.rk = rk(f);

                % System model
                m = kf.models{rk(f)};

                % Predict states at next time instant
                [kf.xkp1_est, kf.Pkp1] = kalman_predict_f(...
                    m.A, m.B, m.Q, ...
                    kf.xk_est, kf.Pk, uk);
                % TODO: For testing only - remove later
                assert(~any(isnan(obj.xkp1_est), 'all'))
                assert(~any(isnan(obj.Pkp1), 'all'))

                % Copy predictions to array
                obj.filters.Xkp1_est(:,:,f) = kf.xkp1_est;
                obj.filters.Pkp1(:,:,f) = kf.Pkp1;

            end
        end
    end
end
