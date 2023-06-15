% Multi-model Kalman Filter class definition
%
% obs = MKFObserverAMM(models,P0,label,x0,p_seq_g_Yk_init,reset)
% Class for simulating a multi-model Kalman filter for state
% estimation of a Markov jump linear system. 
% 
% This version differs from the general MKFObserver in 
% that it models one hypothesis for each system mode with 
% no switching.
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
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
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

classdef MKFObserverAMM < MKFObserver
    methods
        function obj = MKFObserverAMM(models,P0,label,x0, ...
                p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Get number of system models and check their dimensions
            nj = check_models(models);

            % System modes to be modelled
            r0 = (1:nj)';

            % Transition probabilities are all zero
            T = eye(nj);

            % Create super-class observer instance
            obj = obj@MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,false);

            % Store parameters
            obj.type = "MKF_AMM";

            if reset
                % Initialize variables
                obj.reset()
            end

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk, rk)
        % updates the estimates of the multi-model Kalman filter
        % and calculates the predictions of the states and output
        % at the next sample time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %

            % Check size of arguments passed
            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

            % Vector of system modes does not change
            rk = obj.r0;

            % Call update method of super class
            update@MKFObserver(obj, yk, uk, rk);

        end
    end
end
