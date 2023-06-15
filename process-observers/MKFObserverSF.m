% Multi-model Kalman Filter class definition
%
% obs = MKFObserverSF(models,P0,seq,T,label,x0,r0,p_seq_g_Yk_init,reset)
%
% Class for simulating a sub-optimal multi-model observer with 
% sequence fusion, for state estimation of a Markov jump linear 
% system. This version produces posterior estimates of the 
% states and outputs at the current time instant given the data
% available at the current time:
%
%  x_hat(k|k) : estimate of states at time k
%  y_hat(k|k) : estimate of outputs at time k
%
% Prior estimates of the states and outputs at the next
% time instant are also calculated:
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
%   seq : model indicator sequences for each filter (in rows).
%   T : Transition probabity matrix for the Markov switching
%       process. T(i,j) represents the probability of the
%       system switching from model i to j.
%   label : string
%       Arbitrary name to identify the observer.
%   x0 : (n, 1) double
%       Initial state estimates (optional, default zeros).
%   r0 : (1, 1) or (nh, 1) integer (optional)
%       Integer scalar or vector with values in the range 
%       {1, ..., nj} which indicate the prior system modes at time 
%       k = -1. If not provided, the default initialization based 
%       on the mode sequence that is generated will be used.
%   p_seq_g_Yk_init : (nh, 1) double (optional, default uniform)
%       Initial prior probabilities of each hypothesis at
%       time k-1. If not specified, default is equal, i.e.
%       uniform, probability assigned to each hypothesis.
%   reset : logical (default, true)
%       If true, the object's reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       reseting).
%

classdef MKFObserverSF < MKFObserverBM
    properties
        seq (:, 1) cell
        nf (1, 1) double {mustBeInteger}
        i int16
        i_next int16
        idx_branch (1, :) cell
        idx_modes (1, :) cell
        idx_merge (1, :) cell
        nh_max double
    end
    methods
        function obj = MKFObserverSF(models,P0,seq,T,label,x0,r0, ...
            p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                seq (:, 1) cell
                T double
                label (1, 1) string = ""
                x0 = []
                r0 (:, 1) double {mustBeInteger, ...
                    mustBeGreaterThanOrEqual(r0, 1)} = []
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Get number of system models and check their dimensions
            nj = check_models(models);

            % Number of merged hypotheses
            nm = size(seq, 1);

            % Generate the vectors of indices which define the
            % hypothesis branching, mode transitions, and merging 
            % steps of the sequence fusion algorithm
            [idx_branch, idx_modes, idx_merge] = ...
                seq_fusion_indices(cell2mat(seq), nj);

            % Determine maximum number of hypotheses (nh may be
            % time-varying).
            nh_max = max(cellfun(@(idx) size(idx, 1), idx_branch));

            % Determine initial values of discrete states
            nh = size(idx_branch{1}, 1);
            if isempty(r0)
                % Construct rkm1 by branching the merged modes
                % from the last position in the sequence using
                % the current branch index.
                rmkm1 = cellfun(@(s) s(:, end), seq);
                r0 = rmkm1(idx_branch{1});
            elseif isscalar(r0)
                r0 = repmat(r0, nh, 1);
            else
                assert(isequal(size(r0), [nh 1]), "ValueError: size(r0)")
            end

            % Create MKF super-class observer instance
            obj = obj@MKFObserverBM(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,false)

            % Add additional variables used by this observer
            obj.seq = cellfun(@(x) int16(x), seq, 'UniformOutput', false);
            % TODO: allow sequences of different lengths?
            obj.nf = size(cell2mat(seq), 2);
            obj.nm = nm;
            obj.nh_max = nh_max;
            obj.idx_branch = idx_branch;
            obj.idx_modes = idx_modes;
            obj.idx_merge = idx_merge;
            obj.type = "MKF_SF";

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
            reset@MKFObserverBM(obj);

            % Reset sequence index
            obj.i = int16(obj.nf);  % index of last position
            obj.i_next = int16(1);  % index of first position

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

            % Call update method of super class object with the
            % pre-determined branching, mode transitions, and 
            % merging indeces.
            % Note the order of merging, branching, and mode transitions 
            % used here:
            %  - idx_modes(k) are the modes to use for the KF prediction
            %    and update steps at time k to calculate x_est_f(k|k), 
            %    y_est(k|k), and P(k|k), and also to calculate the
            %    transition probabilities.
            %  - idx_merge(k) is then used in the current time step to 
            %    merge  the updated estimates, covariances, using the 
            %    conditoinal probabilities, Pr(Gamma_f(k)|Y(k)).
            %  - Finally, idx_branch(k+1) is used in the current time 
            %    step to branch the merged hypotheses before the KF 
            %    prediction step. However, it must match the modes and 
            %    merging indeces that will occur in the next time step.
            %    Therefore,idx_branch(k+1) is provided here.
            %
            update@MKFObserverBM(obj, yk, uk, ...
                obj.idx_modes{obj.i}, ...
                obj.idx_merge{obj.i}, ...
                obj.idx_branch{obj.i_next});

        end
    end
end
