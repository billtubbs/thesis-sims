% Kalman Filter class definition
%
% obs = KalmanFilterP(model,P0,label,x0,reset)
% Class for simulating a dynamic Kalman filter (i.e.
% with time-varying gain and estimation error
% covariance). This is the prediction form of the KF, 
% which produces predictions of the states and outputs
% at the next time instant given the data at the 
% current time:
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
%       If true, the objects reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       resetting).
%

classdef KalmanFilterP < AbstractLinearFilter
    properties
        P0 double
        K double
        Pkp1 double
    end
    methods
        function obj = KalmanFilterP(model,P0,label,x0,reset)
            arguments
                model struct
                P0 double
                label (1, 1) string = ""
                x0 (:, 1) double = []
                reset logical = true
            end

            % Call super-class constructor
            obj = obj@AbstractLinearFilter(model,"KFP",label,x0,false)

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

            % Note: At initialization time k = 0, the gain is
            % not asigned. It will be calculated at first update.
            obj.K = nan(obj.n, 1);

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk) updates the gain and covariance matrix
        % of the Kalman filter and calculates the estimates of the
        % states and output at the next sample time given the 
        % current measurement and inputs.
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %

            % Update correction gain and covariance matrix
            [obj.K, obj.Pkp1] = kalman_update(obj.Pkp1, obj.model.A, ...
                obj.model.C, obj.model.Q, obj.model.R);

            % Update predictions of states and outputs in 
            % next timestep
            obj.xkp1_est = obj.model.A * obj.xkp1_est ...
                + obj.model.B * uk ...
                + obj.K * (yk - obj.model.C * obj.xkp1_est);
            obj.ykp1_est = obj.model.C * obj.xkp1_est;

        end
    end
end