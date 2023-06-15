% Multi-model Kalman Filter class definition
%
% obs = MKFObserverSP_DI(models,P0,T,d,nh,n_min,label,x0,r0, ...
%     p_seq_g_Yk_init,reset)
%
% Class for simulating a multi-model observer that uses a 
% sequence pruning algorithm for state estimation with a 
% switching system as described in (Andersson, 1985).
% This version includes the detection interval procedure
% described by Robertson et al. (1998) which reduces the
% frequency of the branching and pruning procedure. 
%
% Uses the sequence pruning method described in Eriksson and
% Isaksson (1996):
% - At the updates, let only the most probable sequence,
%   i.e. with the largest weight of all the sequences
%   split into 2 branches.
% - After the update is completed, remove the sequence
%   with the smallest weight.
% - Renormalize the remaining weights to have unit sum.
%
% Restriction to above rules:
% - Do not cut branches immediately after they are born.
%   Let there be a certain minimum life length for all 
%   branches.
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
%   d : integer double
%      Detection interval length in number of sample periods.
%   nh : integer double
%       Number of hypotheses to model (each with a separate
%       Kalman filter).
%   n_min : integer double
%       Minimum life of cloned filters in number of sample
%       periods.
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
% NOTE:
% - The adaptive forgetting component of the AFMM
%   (Andersson, 1985) is not yet implemented.
%
% References:
%  - Eriksson, P.-G., & Isaksson, A. J. (1996). Classification
%    of Infrequent Disturbances. IFAC Proceedings Volumes, 29(1), 
%     6614-6619. https://doi.org/10.1016/S1474-6670(17)58744-3
%  - Andersson, P. (1985). Adaptive forgetting in recursive
%    identification through multiple models. International
%    Journal of Control, 42(5), 1175-1193. 
%    https://doi.org/10.1080/00207178508933420
%  -  Robertson, D. G., & Lee, J. H. (1998). A method for the
%     estimation of infrequent abrupt changes in nonlinear 
%     systems. Automatica, 34(2), 261-270.
%     https://doi.org/10.1016/S0005-1098(97)00192-1%
%

classdef MKFObserverSP_DI < MKFObserverSP
    properties
        d (1, 1) double {mustBeInteger, mustBeNonnegative}
        id int16
        id_next int16
    end
    methods
        function obj = MKFObserverSP_DI(models,P0,T,d,nh,n_min,label, ...
                x0,r0,p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                T double
                d (1, 1) double {mustBeInteger, ...
                    mustBeGreaterThanOrEqual(d,1)}
                nh (1, 1) double {mustBeInteger}
                n_min (1, 1) double {mustBeInteger}
                label (1, 1) string = ""
                x0 = []
                r0 (:, 1) double {mustBeInteger, ...
                    mustBeGreaterThanOrEqual(r0, 1)} = 1
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Create MKF super-class observer instance
            obj = obj@MKFObserverSP(models,P0,T,nh,n_min,label,x0,...
                r0,p_seq_g_Yk_init,false)

            % Add additional variables used by SP observer
            obj.d = d;
            obj.type = "MKF_SP_DI";

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
            reset@MKFObserverSP(obj);

            % Set estimate covariances to high values for all the
            % filters except the first
            for i = 2:obj.nh
                obj.filters.Pkp1(:,:,i) = 1e10 * eye(obj.n);
            end

            % Reset counter for current position in 
            % detection interval.
            obj.id = int16(0);
            obj.id_next = int16(1);

        end
        function update(obj, yk, uk)
        % obs.update(yk, uk)
        % updates the multi-model Kalman filter and calculates the
        % estimates of the states and output at the next sample
        % time.
        %
        % Arguments:
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %

            % Increment sequence index and update counter:
            % obj.id is the counter for the detection interval 
            % (1 <= id <= d). Switching is assumed to only 
            % occur at the end of the detection interval.
            % When obj.id exceeds d (the spacing parameter), it 
            % is reset to 1.
            obj.id = obj.id_next;
            obj.id_next = mod(obj.id, obj.d) + 1;

            % If at start of a detection interval, carry out
            % sequence pruning procedure and update conditional 
            % probability estimates, otherwise, only carry out
            % conditional probability updates.
            %TODO: This needs updating as per MKFObserverSF_DI
            if obj.id == 1

                % Run MKF super-class updates and probability calcs
                update@MKFObserverSP(obj, yk, uk);

            else

                % Run standard MKF probability calcs with same
                % mode transition as at previous time.
                update@MKFObserver(obj, yk, uk, obj.rk, obj.rkm1);

            end

        end
    end
end
