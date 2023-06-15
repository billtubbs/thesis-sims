% Kalman Filter class definition
%
% obs = KalmanFilterF(model,P0,label,x0,reset)
% Class for simulating a dynamic Kalman filter
% (i.e. with time-varying gain and estimation error
% covariance). This is the filtering form of the KF, 
% which produces posterior estimates of the states and
% outputs at the current time instant given the data
% at the current time:
%
%   x_hat(k|k) : estimate of states at time k
%   y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next
% time instant given the data at the current time are
% also calculated:
%
%   x_hat(k+1|k) : estimate of states at time k + 1
%   y_hat(k+1|k) : estimate of outputs at time k + 1
%
% The system model is defined as:
%
%   x(k+1) = A(k) x(k) + B(k) u(k) + w(k)
%     y(k) = C(k) x(k) + v(k)
%
% Note: there is no direct transmission (D = 0).
%
% Arguments:
%   model : struct
%       Struct containing the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and the sampling period, Ts.
%   P0 : matrix, size (n, n)
%       Initial value of covariance matrix of the state
%       estimates at time k = 0.
%   label : string (optional)
%       Name.
%   x0 : vector, size(n, 1), (optional)
%       Initial state estimates at time k = 0.
%   reset : logical (default, true)
%       If true, the object's reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       resetting).
%

classdef KalmanFilterF < AbstractLinearFilter
    properties
        P0 double
        Kf double
        Pkp1 double
        Pk double
        Sk double
    end
    methods
        function obj = KalmanFilterF(model,P0,label,x0,reset)
            arguments
                model struct
                P0 double
                label (1, 1) string = ""
                x0 (:, 1) double = []
                reset logical = true
            end

            % Call super-class constructor
            obj = obj@AbstractLinearFilter(model,"KFF",label,x0,false)

            % Check size of other parameters
            assert(isequal(size(model.Q), [obj.n obj.n]), ...
                "ValueError: size(model.Q)")
            assert(isequal(size(model.R), [obj.ny obj.ny]), ...
                "ValueError: size(model.R)")
            assert(isequal(size(P0), [obj.n obj.n]), ...
                "ValueError: size(P0)")

            % Store parameters
            obj.P0 = P0;

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
            reset@AbstractLinearFilter(obj)

            % Initialize state error covariance P(k|k-1)
            obj.Pkp1 = obj.P0;

            % Note: At initialization time k = 0, these variables
            % are not asigned. They will be calculated at first 
            % update.
            obj.Pk = nan(obj.n);
            obj.Sk = nan(obj.ny);
            obj.Kf = nan(obj.n, obj.ny);

        end
        function correct(obj, yk)
        % Calculate updated state estimates x_est(k|k) output 
        % estimates y_est(k|k) and error covariance P(k|k) at 
        % the current time using the current output measurement 
        % y(k). The correction gain Kf(k) and the covariance of 
        % the output estimation error S(k) are also calculated.

            [obj.xk_est, obj.Pk, obj.yk_est, obj.Kf, obj.Sk] = ...
                kalman_update_f(obj.model.C, obj.model.R, ...
                obj.xkp1_est, obj.Pkp1, yk);

        end
        function predict(obj, uk)
        % Calculate predicted states x_est(k+1|k) and outputs, 
        % y(k+1|k) at next time instant using the system model and 
        % the current control input, u(k).

            [obj.xkp1_est, obj.Pkp1] = kalman_predict_f(obj.model.A, ...
                obj.model.B, obj.model.Q, obj.xk_est, obj.Pk, uk);
            obj.ykp1_est = obj.model.C * obj.xkp1_est;

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) carries out the following steps:
        %  1. Calculates the correction gain, Kf(k) and 
        %     output estimation error covariance, S(k).
        %  2. Corrects the estimates of the states and outputs at
        %     the current time using the current measurement, y(k):
        %         x_est(k|k-1) -> x_est(k|k)
        %         y_est(k|k-1) -> y_est(k|k)
        %  3. Calculates the predicted states and outputs at the
        %     next time instant:
        %         x_est(k+1|k)
        %         y_est(k+1|k)
        %     given the current input, u(k), and the system model.
        %
        % To do steps 1 and 2 only use: obs.correct(yk).
        % To do step 3 only use: obs.predict(uk).
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %
        % After calling this method, the following properties will
        % be updated:
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

            % Update estimates based on current measurement
            obj.correct(yk)

            % Predictions based on current inputs
            obj.predict(uk)

        end
    end
end