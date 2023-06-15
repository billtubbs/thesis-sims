% Multi-model Kalman Filter class definition
%
% obs = MKFObserverS(models,P0,seq,T,label,x0, ...
%     p_seq_g_Yk_init,reset)
%
% Class for simulating a multi-model Kalman filter for state
% estimation of a Markov jump linear system. 
% 
% This version differs from MKFObserver in that it
% has a pre-determined set of sequences of system modes
% (i.e. transitions) that it models.
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
%   seq : (nh, 1) cell
%       Set of sequences of system models representing 
%       alternative hypotheses of the mode transitions over
%       a fixed period of time. Each sequence should be a
%       row vector of integers in the range {1, ..., nj} 
%       which indicate the system mode at time k = 0, 1, 
%       ..., nf-1.
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
%       resetting).
%

classdef MKFObserverS < MKFObserver
    properties
        seq (:, 1) cell
        nf (1, 1) double {mustBeInteger}
        i int16
        i_next int16
    end
    methods
        function obj = MKFObserverS(models,P0,seq,T,label,x0, ...
                r0,p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                seq (:, 1) cell
                T double
                label (1, 1) string = ""
                x0 = []
                r0 (:, 1) int16 {mustBeGreaterThan(r0, 0)} = 1
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Determine initial system modes at time k = -1
            if isscalar(r0)
                r0 = repmat(r0, size(seq, 1), 1);
            else
                assert(isequal(size(r0), [size(seq, 1) 1]), ...
                    "ValueError: r0")
            end

            % Create super-class observer instance
            obj = obj@MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,false);

            % Add additional variables used by SP observer
            obj.seq = cellfun(@(x) int16(x), seq, 'UniformOutput', false);
            obj.nf = size(cell2mat(seq), 2);  % TODO: allow sequences of
                                              %     of different lengths
            obj.type = "MKF_S";

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

            % Call reset method of super class object
            reset@MKFObserver(obj);

            % Reset sequence index
            obj.i = int16(0);
            obj.i_next = int16(1);

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
        % updates the estimates of the multi-model Kalman filter
        % and calculates the predictions of the states and output
        % at the next sample time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %

            % Increment sequence index (at end of sequence it 
            % wraps to beginnning)
            obj.i = obj.i_next;
            obj.i_next = mod(obj.i, obj.nf) + 1;

            % Get vector of current system modes from sequence
            obj.rk = cellfun(@(s) s(:, obj.i), obj.seq);

            % Call update method of super class object
            update@MKFObserver(obj, yk, uk, obj.rk);

        end
    end
end
