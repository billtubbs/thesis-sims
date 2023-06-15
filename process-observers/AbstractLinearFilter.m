% Abstract class definition for linear filters
%
% obs = AbstractLinearFilter(model,type,label,x0,reset)
%
% Abstract Class for inheriting when defining other
% filter classes which use a linear state-space model, 
% such as the Kalman and Luenberger filters.
%
% Arguments:
%   model : struct
%       Struct containing the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and the sampling period, Ts.
%   type : string
%       Code used to categorize the specific type of 
%       filter.
%   label : string (optional)
%       Name.
%   x0 : vector, size(n, 1), (optional)
%       Initial state estimates at time k = 0.
%   reset : logical (default, true)
%       If true, the object's reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       reseting).
%

classdef (Abstract) AbstractLinearFilter < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        model struct
        Ts (1, 1) double {mustBeNonnegative}
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        xk_est (:, 1) double
        yk_est (:, 1) double
        label (1, 1) string
        x0 (:, 1) double
        type (1, 1) string  % is this still needed? use classdef
    end
    methods
        function obj = AbstractLinearFilter(model,type,label,x0,reset)
            arguments
                model struct
                type (1, 1) string
                label (1, 1) string = ""
                x0 (:, 1) double = []
                reset logical = true
            end

            % Check model struct is specified correctly
            [n, nu, ny, Ts, ~] = check_model(model);

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial values
            else
                assert(isequal(size(x0), [n 1]), ...
                    "ValueError: size(x0)")
            end
            if label == ""
                label = type;
            end

            % Store parameters
            obj.n = n;
            obj.nu = nu;
            obj.ny = ny;
            obj.model = model;
            obj.Ts = Ts;
            obj.label = label;
            obj.x0 = x0;
            obj.type = type;

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
            % Note: At initialization time k = 0, xkp1_est and
            % ykp1_est represent prior estimates of the states,
            % i.e. x_est(k|k-1) and y_est(k|k-1).
            obj.xkp1_est = obj.x0;
            obj.ykp1_est = obj.model.C * obj.xkp1_est;

            % At initialization time k = 0, x_est(k|k)
            % and y_est(k|k) have not yet been computed.
            obj.xk_est = nan(obj.n, 1);
            obj.yk_est = nan(obj.ny, 1);

        end
    end
end