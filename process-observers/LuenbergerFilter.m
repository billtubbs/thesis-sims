% Luenberger filter class definition
%
% obs = LuenbergerFilter(model,poles,label,x0,reset)
% Class for simulating a steady-state Kalman filter
% (i.e. with static gain).
%
% Arguments:
%   model : struct
%       Struct containing the parameters of a linear
%       model of the system dynamics. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and the sampling period, Ts.
%   poles : vector
%       Desired closed-loop poles of filter dynamics.
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
% References:
%  -  D. Luenberger, "An introduction to observers," in IEEE 
%     Transactions on Automatic Control, vol. 16, no. 6, 
%     pp. 596-602, December 1971, doi: 10.1109/TAC.1971.1099826.
%

classdef LuenbergerFilter < AbstractLinearFilter
    properties
        poles {mustBeNumeric}
        K double
        Pkp1 double
    end
    methods
        function obj = LuenbergerFilter(model,poles,label,x0,reset)
            arguments
                model struct
                poles (:, 1) double
                label (1, 1) string = ""
                x0 (:, 1) double = []
                reset logical = true
            end

            % Call super-class constructor
            obj = obj@AbstractLinearFilter(model,"LB",label,x0,false)

            % Set properties for Luenberger filter
            obj.poles = poles;

            % Compute observer gain
            if obj.ny == 1
                obj.K = acker(obj.model.A', obj.model.C', poles)';
            else
                obj.K = place(obj.model.A', obj.model.C', poles)';
            end

            if reset
                % Initialize variables
                obj.reset()
            end

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk) updates the estimates of the
        % states and output at the next sample time:
        %         x_est(k+1|k)
        %         y_est(k+1|k)
        % given the current measurement y(k) and the current 
        % input u(k).
        %
        % Arguments:
        %   yk : vector, size (ny, 1)
        %       System output measurements at current time k.
        %   uk : vector, size (nu, 1), optional
        %       System inputs at current time k.
        %
        % After calling this method, the following properties will
        % be updated:
        %   - obs.xkp1_est  (x_est(k+1|k))
        %   - obs.ykp1_est  (y_est(k+1|k))
        %

            % Update predictions of states and outputs in 
            % next timestep
            obj.xkp1_est = obj.model.A * obj.xkp1_est ...
                + obj.model.B * uk ...
                + obj.K * (yk - obj.model.C * obj.xkp1_est);
            obj.ykp1_est = obj.model.C * obj.xkp1_est;

        end
    end
end
