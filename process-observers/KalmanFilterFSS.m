% Steady-state Kalman Filter class definition
%
% obs = KalmanFilterFSS(model,label,x0,reset)
% Class for simulating a steady-state Kalman filter
% (i.e. with static gain). This is the filtering form of 
% the KF, which produces posterior estimates of the 
% states and outputs at the current time instant given 
% the data at the current time:
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
%   label : string (optional)
%       Name.
%   x0 : vector, size(n, 1), (optional)
%       Initial state estimates at time k = 0.
%   reset : logical (default, true)
%       If true, the objects reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       resetting).
%

classdef KalmanFilterFSS < AbstractLinearFilter
    properties
        P0 double
        Kf double
        Pkp1 double
        Sk double
    end
    methods
        function obj = KalmanFilterFSS(model,label,x0,reset)
            arguments
                model struct
                label (1, 1) string = ""
                x0 (:, 1) double = []
                reset logical = true
            end

            % Call super-class constructor
            obj = obj@AbstractLinearFilter(model,"KFFSS",label,x0,reset)

            % Check size of other parameters
            assert(isequal(size(model.Q), [obj.n obj.n]), ...
                "ValueError: size(model.Q)")
            assert(isequal(size(model.R), [obj.ny obj.ny]), ...
                "ValueError: size(model.R)")

            % Compute the steady-state gain and error covariance matrix
            % This is the gain for the filtering form of the KF:
            [obj.Kf, obj.Pkp1] = kalman_gain_ss(obj.model.A, obj.model.C, ...
                obj.model.Q, obj.model.R);

            % Error covariance of output estimation error
            obj.Sk = obj.model.C * obj.Pkp1 * obj.model.C' + obj.model.R;

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

            % At initialization time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);

        end
        function correct(obj, yk)
        % Calculate updated state estimates x_est(k|k) output 
        % estimates y_est(k|k) at the current time using the 
        % current output measurement y(k).

            % Update prior state estimates using measurements 
            % from current time to produce 'a posteriori' state
            % estimates
            obj.xk_est = obj.xkp1_est + obj.Kf * (yk - obj.ykp1_est);

            % Updated output estimate
            obj.yk_est = obj.model.C * obj.xk_est;

        end
        function predict(obj, uk)
        % Calculate predicted states x_est(k+1|k) and outputs, 
        % y(k+1|k) at next time instant using the system model and 
        % the current control input, u(k).

            obj.xkp1_est = obj.model.A * obj.xk_est + obj.model.B * uk;
            obj.ykp1_est = obj.model.C * obj.xkp1_est;

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) carries out the following steps:
        %  1. Corrects the estimates of the states and outputs at
        %     the current time using the current measurement, y(k):
        %         x_est(k|k-1) -> x_est(k|k)
        %         y_est(k|k-1) -> y_est(k|k)
        %  3. Calculates the predicted states and outputs at the
        %     next time instant:
        %         x_est(k+1|k)
        %         y_est(k+1|k)
        %     given the current input u(k) and the system model.
        %
        % To do step 1 only use: obs.correct(yk).
        % To do step 2 only use: obs.predict(yk).
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %
        % After calling this method, the following properties will
        % be updated:
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