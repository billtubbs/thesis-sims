% Switching Kalman Filter class definition
%
% obs = SKFObserver(models,P0,label,x0,r0,reset)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance) for state estimation of a jump linear
% system where the system model switches between
% a finite set of models each time step.
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
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
%   r0 : (nh, 1) integer (optional, default ones)
%       Integer in the range {1, ..., nj} which indicates
%       the prior system mode at time k = -1.
%   reset : logical (default, true)
%       If true, the objects reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       resetting).
%

classdef SKFObserver < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        nu (1, 1) double {mustBeInteger}
        ny (1, 1) double {mustBeInteger}
        n (1, 1) double {mustBeInteger}
        nj (1, 1) double {mustBeInteger}
    end
    properties
        models (1, :) cell
        Ts (1, 1) double {mustBeNonnegative}
        Kf double
        P0 double
        label (1, 1) string
        x0 (:, 1) double
        r0 (1, 1) int16
        xk_est (:, 1) double
        Pk double
        yk_est (:, 1) double
        xkp1_est (:, 1) double
        Pkp1 double
        Sk double
        rk (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(rk, 1)}
        rkm1 (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(rkm1, 1)}
        type (1, 1) string  % is this still needed? use classdef
    end
    methods
        function obj = SKFObserver(models,P0,label,x0,r0,reset)
            arguments
                models (1, :) cell
                P0 double
                label (1, 1) string = ""
                x0 (:, 1) double = []
                r0 (1, 1) {mustBeInteger} = 1
                reset logical = true
            end

            % Get number of system models and check their dimensions
            [nj, n, nu, ny, Ts] = check_models(models);

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
                assert(isequal(size(x0), [n 1]))
            end

            % Store parameters
            obj.Ts = Ts;
            obj.nj = nj;
            obj.nu = nu;
            obj.n = n;
            obj.ny = ny;
            obj.models = models;
            obj.P0 = P0;
            obj.x0 = x0;
            obj.r0 = r0;
            obj.type = "SKF";
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

            % Initialize state estimate and system mode
            % Note: At initialization at time k = 0, xkp1_est and
            % r0 represent prior estimates of the states and mode,
            % i.e. x_est(k|k-1) and r(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.rk = obj.r0;
            obj.rkm1 = [];

            % Initialize estimate covariance
            obj.Pkp1 = obj.P0;

            % At initialization time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);
            obj.Pk = nan(obj.n);
            obj.Sk = nan(obj.ny);

            % Gain will be calculated at first update
            obj.Kf = nan(obj.n, obj.ny);

        end
        function update_step(obj, yk, rk)
        % Calculate updated state estimates, x_est(k|k), and error 
        % covariance, P(k|k) at the current time using the current
        % output measurement, y(k) and the current system model, 
        % r(k). The correction gain, Kf(k) and covariance of the 
        % output estimation error, S(k) are also calculated.
            obj.rk = rk;
            [obj.xk_est, obj.Pk, obj.yk_est, obj.Sk, obj.Kf] = ...
                kalman_update_f(obj.models{rk}.C, ...
                obj.models{rk}.R, obj.xkp1_est, obj.Pkp1, yk);
        end
        function prediction_step(obj, uk, rk)
        % Calculate predicted states, x_est(k+1|k), and outputs, 
        % y(k+1|k) at next time instant using the current system 
        % model, r(k), and the current control input, u(k).
            obj.rk = rk;
            [obj.xkp1_est, obj.Pkp1] = kalman_predict_f(...
                obj.models{rk}.A, obj.models{rk}.B, obj.models{rk}.Q, ...
                obj.xk_est, obj.Pk, uk);
        end
        function update(obj, yk, uk, rk)
        % obs.update(yk, uk) carries out the following steps:
        %  1. Calculates the correction gain, Kf(k) and 
        %     output estimation error covariance, S(k).
        %  2. Updates the estimates of the states and outputs at
        %     the current time using the current measurement, y(k):
        %         x_est(k|k-1) -> x_est(k|k)
        %         y_est(k|k-1) -> y_est(k|k)
        %  3. Calculates the predicted states and outputs at the
        %     next time instant:
        %         x_est(k+1|k)
        %         y_est(k+1|k)
        %     given the current input, u(k), and the current system
        %     model.
        %
        % To do steps 1 and 2 only use: obs.update_step(yk, rk).
        % To do step 3 only use: obs.prediction_step(yk, rk).
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %
        % After calling this method, the following properties will
        % be updated:
        %   - obs.rk  (r(k))
        %   - obs.Kf  (Kf(k))
        %   - obs.Sk  (S(k))
        %   - obs.xk_est  (x_est(k|k))
        %   - obs.yk_est  (y_est(k|k))
        %   - obs.xkp1_est  (x_est(k+1|k))
        %   - obs.ykp1_est  (y_est(k+1|k))
        %

            % Check size of arguments passed
            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")
            assert(isequal(size(rk), [1 1]), "ValueError: size(rk)")

            % Update estimates based on current measurement 
            % and system mode
            obj.update_step(yk, rk)

            % Predictions based on current inputs and system 
            % mode
            obj.prediction_step(uk, rk)

        end
    end
end