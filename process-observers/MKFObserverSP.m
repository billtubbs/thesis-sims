% Multi-model Kalman Filter class definition
%
% obs = MKFObserverSP(models,P0,T,nh,n_min,label,x0,r0, ...
%     p_seq_g_Yk_init,reset)
%
% Object class for simulating a multi-model observer that
% uses a sequence pruning algorithm for state estimation 
% with a switching system as described in (Andersson, 1985).
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
%       If true, the object's reset method is called after
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
%

classdef MKFObserverSP < MKFObserver
    properties (SetAccess = immutable)
        n_min double {mustBeInteger, mustBeNonnegative}
        n_main double
        n_hold double
    end
    properties
        f_main
        f_hold
    end
    methods
        function obj = MKFObserverSP(models,P0,T,nh,n_min,label,x0,r0, ...
            p_seq_g_Yk_init,reset)
            arguments
                models (1, :) cell
                P0 double
                T double
                nh (1, 1) double {mustBeInteger, ...
                    mustBeGreaterThan(nh, 0)}
                n_min (1, 1) double {mustBeInteger, ...
                    mustBeGreaterThan(n_min, 0)}
                label (1, 1) string = ""
                x0 = []
                r0 (:, 1) double {mustBeInteger, ...
                    mustBeGreaterThan(r0, 0)} = 1
                p_seq_g_Yk_init = []
                reset logical = true
            end

            % Number of system modes and model dimensions
            [nj, n, ~, ~] = check_models(models);

            % Determine initial values of discrete states
            % TODO: Is r0 needed for this observer?
            if isscalar(r0)
                % By default, assume initial hypothesis sequences 
                % at time k = -1 all start with 1 (i.e. no shocks)
                r0 = repmat(r0, nh, 1);
            else
                assert(isequal(size(r0), [nh 1]), "ValueError: r0")
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]), "ValueError: x0")
            end

            % Define filter groups ('main', 'holding' and 'unused')
            n_min = int16(n_min);
            assert(n_min > 0)
            assert(nh > 0)
            n_hold = (nj - 1) * n_min;
            n_main = nh - n_hold;

            % Check there are enough filters in total to accommodate
            % those in the holding group + at least one in main group
            assert(n_main >= nj - 1, "ValueError: nh is too low.")

            % Create MKF super-class observer instance
            obj = obj@MKFObserver(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,false);

            % Add additional variables used by this observer
            obj.n_min = n_min;
            obj.n_main = n_main;
            obj.n_hold = n_hold;
            obj.type = "MKF_SP";

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

            % Set estimate covariances to high values for all the
            % filters except the first
            % TODO: Do we need this for SP algorithm?
            %for i = 2:obj.nh
            %    obj.filters.Pkp1(:,:,i) = 1e10 * eye(obj.n);
            %end

            % Set initial filter indices
            obj.f_main = int16(1:obj.n_main);
            obj.f_hold = int16(obj.n_main+1:obj.n_main+obj.n_hold);

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

            % Implementation of filter pruning algorithm
            %
            % - As described in Eriksson and Isaksson (1996):
            %
            %   "Do not cut branches immediately after they are 
            %    born. Let there be a certain minimum life length 
            %    for all branches."
            %
            % To implement this rule, filters are organized into
            % two groups:
            % 
            %  1. Main group : obj.filters{f} for f in obj.f_hold
            %     Longer-surviving filters which can be killed off.
            %  2. Holding group : obj.filters{f} for f in obj.f_main
            %     Recently created filters which cannot be killed 
            %     off and are held for n_min time periods before
            %     being transferred to the main group.
            %
            % Note that for n_dist input disturbances, n_dist new
            % sequences are created each time step, so the
            % holding group needs to be n_dist*obj.n_min in size.

            % Index of current most likely sequence
            [~, f_max] = max(obj.p_seq_g_Yk);

            % Number of new branches to add
            n_new = obj.nj - 1;

            % Left-shift all filters in holding group. This causes
            % the n_new left-most values to 'roll-over' to the right
            % of f_hold. E.g. 
            % circshift([1 2 3], -1) -> [2 3 1]
            obj.f_hold = circshift(obj.f_hold, -n_new);

            % Filters to be moved out of holding group
            f_move = obj.f_hold(end-n_new+1:end);

            % Rank sequences in main group according to 
            % conditional probabilities
            [~, i_rank] = sort(obj.p_seq_g_Yk(obj.f_main), 'descend');
            obj.f_main = obj.f_main(i_rank);

            % Select those with lowest probability for pruning
            f_to_prune = obj.f_main(end-n_new+1:end);

            % Replace pruned sequences with those from holding
            % group
            obj.f_main(end-n_new+1:end) = f_move;

            % Set all mode indicator values to 1 (no shock)
            rk = ones(obj.nh, 1);

            % Make clones of most probable hypotheses and filters
            % and put new filters at start of holding group
            obj.f_hold(end-n_new+1:end) = f_to_prune;
            for i = 1:n_new
                obj.copy_filters(f_max, f_to_prune(i))
                % Change mode indicator of cloned filter to 
                % shock-occurence mode
                rk(f_to_prune(i)) = i + 1;
            end

            % TODO: Add online noise variance estimation with
            % forgetting factor as described in Anderson (1995)
            % equation (4.3).

            % Run MKF super-class updates and probability calcs
            update@MKFObserver(obj, yk, uk, rk);

        end
        function copy_filters(obj, a, b)
        % Make a copy of filter in position a and save it
        % in position b.

            % Transfer filter variables from a to b
            obj.filters.Xkp1_est(:, :, b) = obj.filters.Xkp1_est(:, :, a);
            obj.filters.Pkp1(:, :, b) = obj.filters.Pkp1(:, :, a);
            obj.filters.Xk_est(:, :, b) = obj.filters.Xk_est(:, :, a);
            obj.filters.Pk(:, :, b) = obj.filters.Pk(:, :, a);
            obj.filters.Yk_est(:, :, b) = obj.filters.Yk_est(:, :, a);

            % Transfer probability information
            obj.p_seq_g_Yk(b) = obj.p_seq_g_Yk(a);
            obj.p_rk_g_rkm1(b) = obj.p_rk_g_rkm1(a);
            obj.p_seq_g_Ykm1(b) = obj.p_seq_g_Ykm1(a);

            % Switch sequence histories - not yet implemented
            % TODO: Implement in parent class MKF_SP_RODD
            %obj.seq{f_to_prune(i)} = obj.seq{f_max};
            % Set next sequence value to index of the shock
            %obj.seq{f_to_prune(i)}(:, obj.i_next(1)) = i + 1;

        end
    end
end
