% Multi-model Kalman Filter class definition
%
% Object class for simulating a multi-model observer
% using the adaptive forgetting through multiple models 
% (AFMM) algorithm for state estimation in the presence of 
% infrequently-occurring deterministic disturbances, as 
% described in Eriksson and Isaksson (1996). This version
% differs from MKFObserverSP in that it uses instances of
% KalmanFilterF to do the hypothesis filter updates and
% predictions.
%
% Uses a sequence pruning method described in Eriksson and
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
% obs = MKFObserverSPF(model,u_meas,P0,epsilon, ...
%     sigma_wp,Q0,R,nh,n_min,label,x0,r0,reset)
%
% Creates a struct for simulating a multi-model observer
% using the adaptive forgetting through multiple models 
% (AFMM) method for state estimation in the presence of 
% infrequently-occurring deterministic disturbances, as 
% described in Eriksson and Isaksson (1996).
%
% Arguments:
%   model : struct
%       Struct containing the parameters of a linear
%       model of the system dynamics including disturbances
%       and unmeasured inputs. These include: A, B, 
%       and C for the system matrices, and the sampling 
%       period, Ts.
%   u_meas : (:, 1) logical
%       Binary vector indicating which inputs are measured.
%   P0 : (n, n) double
%       Initial value of covariance matrix of the state
%       estimates.
%   epsilon : (nw, 1) double
%       Probability(s) of shock disturbance(s).
%   sigma_wp : (nw, 2) double
%       Standard deviation(s) of shock disturbance(s).
%   Q0 : (n, n) double
%       Process noise covariance matrix (n, n) with 
%       variances for each state on the diagonal. The  
%       values for states impacted by the unmeasured input
%       disturbances should be set to zero as the
%       appropriate variances will be added by the
%       algorithm during observer updates.
%   R : (ny, ny) double
%       Output measurement noise covariance matrix.
%   nh : integer double
%       Number of models (Kalman filters) to utilise.
%   n_min : integer double
%       Minimum life of cloned filters in number of sample
%       periods.
%   label : String name.
%   x0 : Initial state estimates (optional, default zeros)
%   reset : logical (default, true)
%       If true, the objects reset method is called after
%       initialization (this is mainly intended for use by
%       other objects instantiating an instance without
%       reseting).
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

classdef MKFObserverSPF < MKFObserverF
    properties (SetAccess = immutable)
        u_meas (:, 1) logical {mustBeNumericOrLogical}
        n_min double {mustBeInteger, mustBeNonnegative}
    end
    properties
        sys_model struct
        %seq (:, 1) cell
        %i (1, 1) {mustBeInteger, mustBeNonnegative}
        %i_next (1, 1) {mustBeInteger, mustBeNonnegative}
        p_seq double
        Q0 double
        R double
        epsilon double
        sigma_wp double
        n_main
        n_hold
        f_main
        f_hold
    end
    methods
        function obj = MKFObserverSPF(model,u_meas,P0,epsilon, ...
                sigma_wp,Q0,R,nh,n_min,label,x0,r0,reset)
            arguments
                model struct
                u_meas (:, 1) logical
                P0 double
                epsilon double
                sigma_wp double
                Q0 (:, :) double
                R (:, :) double
                nh (1, 1) double {mustBeInteger}
                n_min (1, 1) double {mustBeInteger}
                label (1, 1) string = ""
                x0 = []
                r0 (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(r0, 1)} = []
                reset logical = true
            end

            % Number of states
            [n, nu, ny] = check_dimensions(model.A, model.B, model.C);
            if isprop(model, "D")
                assert(isequal(model.D, zeros(ny, nu)), ...
                    "ValueError: direct transmission (model.D)")
            end

            % Check size of initial process covariance matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Determine initial state values
            if isempty(r0)
                % By default, assume initial hypothesis sequences 
                % at time k = -1 all start with 1 (i.e. no shocks)
                r0 = ones(nh, 1);
            else
                assert(isequal(size(r0), [nh 1]), "ValueError: r0")
            end

            % Determine initial state values
            if isempty(x0)
                x0 = zeros(n, 1);  % default initial states
            else
                assert(isequal(size(x0), [n 1]), "ValueError: x0")
            end

            % Observer model without disturbance noise input
            nw = sum(~u_meas);  % Number of input disturbances
            assert(nw > 0, "ValueError: u_meas");
            Bu = model.B(:, u_meas);
            Bw = model.B(:, ~u_meas);

            % Construct process noise covariance matrices for each 
            % possible input disturbance (returns a cell array)
            [Q, p_rk_g_rkm1] = construct_Q_model_SP(Q0, Bw, epsilon, ...
                sigma_wp.^2, nw);

            % Number of models (each with a different hypothesis sequence)
            nj = numel(Q);

            % Transition probability matrix
            % Note that for RODD disturbances Pr(r(k)) is
            % assumed to be an independent random variable.
            T = repmat(p_rk_g_rkm1', nj, 1);

            % Initialize indicator sequences with 1s
            %seq = mat2cell(int16(ones(nh, nf)), int16(ones(1, nh)), nf);

            % Define filter groups ('main', 'holding' and 'unused')
            n_min = int16(n_min);
            assert(n_min > 0)
            assert(nh > 0)
            n_hold = nw * n_min;
            n_main = nh - n_hold;

            % Check there are enough filters in total to accommodate
            % those in the holding group + at least one in main group
            assert(n_main >= nw, "ValueError: nh is too low.")

            % Filter indices
            f_main = int16(1:n_main);
            f_hold = int16(n_main+1:n_main+n_hold);

            % System models are all the same - only Q switches
            model_obs = model;
            model_obs.B = Bu;
            models = repmat({model_obs}, 1, nj);
            for i = 1:nj
                models{i}.Q = Q{i};
                models{i}.R = R;
            end

            % Sequence pruning algorithm initialization
            % Assign all probability to first filter.
            p_seq_g_Yk_init = [1; zeros(nh-1, 1)];

            % Create MKF super-class observer instance
            obj = obj@MKFObserverF(models,P0,T,r0,label,x0, ...
                p_seq_g_Yk_init,false);

            % Add additional variables used by AFMM observer
            obj.sys_model = model;
            obj.u_meas = u_meas;
            %obj.seq = seq;
            obj.n_min = n_min;
            obj.n_main = n_main;
            obj.n_hold = n_hold;
            obj.f_main = f_main;
            obj.f_hold = f_hold;
            obj.P0 = P0;
            obj.Q0 = Q0;
            obj.R = R;
            obj.epsilon = epsilon;
            obj.sigma_wp = sigma_wp;
            obj.type = "MKF_SPF";

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
            reset@MKFObserverF(obj);

            % Set estimate covariances to high values for the
            % rest of the filters
            for i = 2:obj.nh
                obj.filters.objects{i}.Pkp1 = 1e10 * eye(obj.n);
            end

            % Reset sequence index
            %obj.i = int16(0);
            %obj.i_next = int16(1);

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

            % Set next sequence value to 1 (no shock) for all
            % sequences
            %for f = 1:obj.nh
            %    obj.seq{f}(:, obj.i_next(1)) = 1;
            %end

            % Number of disturbances
            nw = size(obj.epsilon, 1);

            % Left-shift all filters in holding group. This causes
            % the last nw values to 'roll-over' to the left of f_hold.
            % e.g. circshift([1 2 3], 1) -> [3 1 2]
            obj.f_hold = circshift(obj.f_hold, -nw);

            % Filters to be moved out of holding group
            f_move = obj.f_hold(end-nw+1:end);

            % Rank sequences in main group according to 
            % conditional probabilities
            [~, i_rank] = sort(obj.p_seq_g_Yk(obj.f_main), 'descend');
            obj.f_main = obj.f_main(i_rank);

            % Select those with lowest probability for pruning
            f_to_prune = obj.f_main(end-nw+1:end);

            % Replace pruned sequences with those from holding
            % group
            obj.f_main(end-nw+1:end) = f_move;

            % Set all mode indicator values to 1 (no shock)
            rk = ones(obj.nh, 1);

            % Make clone(s) of most probable sequence and filter
            % and put new filter(s) at start of holding group
            obj.f_hold(end-nw+1:end) = f_to_prune;
            for i = 1:nw
                label = obj.filters.objects{f_to_prune(i)}.label;
                obj.filters.objects{f_to_prune(i)} = obj.filters.objects{f_max}.copy();
                obj.filters.objects{f_to_prune(i)}.label = label;  % keep label
                obj.p_seq_g_Yk(f_to_prune(i)) = obj.p_seq_g_Yk(f_max);
                obj.p_rk_g_rkm1(f_to_prune(i)) = obj.p_rk_g_rkm1(f_max);
                obj.p_seq_g_Ykm1(f_to_prune(i)) = obj.p_seq_g_Ykm1(f_max);
                %obj.seq{f_to_prune(i)} = obj.seq{f_max};
                % Set next sequence value to index of the shock
                %obj.seq{f_to_prune(i)}(:, obj.i_next(1)) = i + 1;
                % Set mode indicator of cloned filter to shock 
                % occurs value
                rk(f_to_prune(i)) = i + 1;
            end

            % Copy filter variables into arrays
            for f = 1:obj.nh
                obj.filters.Xk_est(:,:,f) = obj.filters.objects{f}.xk_est;
                obj.filters.Pk(:,:,f) = obj.filters.objects{f}.Pk;
                obj.filters.Yk_est(:,:,f) = obj.filters.objects{f}.yk_est;
                obj.filters.Xkp1_est(:,:,f) = obj.filters.objects{f}.xkp1_est;
                obj.filters.Pkp1(:,:,f) = obj.filters.objects{f}.Pkp1;
            end

            % TODO: Add online noise variance estimation with
            % forgetting factor as described in Anderson (1995)
            % equation (4.3).

            % Run MKF super-class updates and probability calcs
            update@MKFObserverF(obj, yk, uk, rk);

        end
    end
end
