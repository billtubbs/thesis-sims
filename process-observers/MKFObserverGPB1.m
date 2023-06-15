% Multi-model Kalman Filter class definition
%
% obs = MKFObserverGPB1(models,P0,T,label,x0,p_seq_g_Yk_init)
%
% Class for simulating the generalised pseudo-Bayes multi-
% model Kalman filter for state estimation of Markov jump
% linear systems. This is the first-order version of the
% algorithm (GPB1).
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
%   T : Transition probabity matrix of the Markov switching
%       process.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
%   p_seq_g_Yk_init : (optional, default uniform)
%       Initial prior hypothesis probabilities at time k-1.
%       If not specified, default is equal, i.e. uniform,
%       probability assigned to each hypothesis.
%

classdef MKFObserverGPB1 < MKFObserver
    methods
        function obj = MKFObserverGPB1(models,P0,T,label,x0, ...
                p_seq_g_Yk_init)
            arguments
                models (1, :) cell
                P0 double
                T double
                label (1, 1) string = ""
                x0 = []
                p_seq_g_Yk_init = []
            end

            % Get number of system models and check their dimensions
            nj = check_models(models);

            % System modes to be modelled
            r0 = (1:nj)';

            % Create super-class observer instance
            obj = obj@MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,false);

            % Store additional parameters
            obj.type = "MKF_GPB1";

            % Initialize all variables
            obj.reset()

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
        % updates the estimates of the switching Kalman filter
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

            % Vector of system modes (same at all timesteps)
            obj.rk = obj.r0;

            % Update state and output estimates based on current
            % measurement and prior predictions
            [obj.xk_est, obj.yk_est, obj.Pk, obj.p_seq_g_Yk] = ...
                GPB1_update( ...
                    obj.models, ...
                    obj.T, ...
                    obj.filters.Xkp1_est, ...
                    obj.filters.Pkp1, ...
                    yk, ...
                    obj.p_seq_g_Yk ...
                );

            % Calculate predictions of each filter in next time instant
            % GBP1 branches the estimates from previous time
            % instant when making predictions for next:
            %   xi_est(k+1|k) = Ai(k) * x_est(k|k-1) + Bi(k) * u(k);
            %   Pi(k+1|k) = Ai(k) * P(k|k-1) * Ai(k)' + Qi(k);
            %
            for j = 1:obj.nh
                m = obj.models{obj.rk(j)};
                [obj.filters.Xkp1_est(:,:,j), obj.filters.Pkp1(:,:,j)] = ...
                    kalman_predict_f(m.A, m.B, m.Q, obj.xk_est, obj.Pk, uk);
            end

        end
    end
end