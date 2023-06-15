% Multi-model Kalman Filter class definition
%
% obs = MKFObserverSP_RODD(model,io,P0,epsilon, ...
%     sigma_wp,Q0,R,nh,n_min,label,x0,r0)
%
% Object class for simulating a multi-model observer that
% uses the adaptive forgetting through multiple models 
% (AFMM) algorithm for state estimation in the presence of 
% infrequently-occurring deterministic disturbances, as 
% described in Eriksson and Isaksson (1996).
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
% The observer object can be used recursively in an 
% iteration loop or in a Simulink S-function block (see 
% MKFObserver_sfunc.m)
%
% Arguments:
%   model : struct
%       Struct containing the parameters of a linear
%       model of the system dynamics including disturbances
%       and unmeasured inputs. These include: A, B, 
%       and C for the system matrices, and the sampling 
%       period, Ts.
%   io : struct
%       Struct containing logical vectors u_known and y_meas
%       indicating which inputs are known/unknown and which 
%       outputs are measured/unmeasured.
%   P0 : (n, n) double
%       Initial value of covariance matrix of the state
%       estimates.
%   epsilon : (nw, 1) double
%       Probability(s) of shock disturbance(s).
%   sigma_wp : (1, nw) cell array
%       Standard deviations of disturbances. Each element of
%       the cell array is either a scalar for a standard (Gaussian)
%       noise or a (1, 2) vector for a random shock disturbance.
%   Q0 : (n, n)
%       Matrix containing variance values for process
%       states. Only values in the rows and columns corresponding 
%       to the process states are used, usually the upper left
%       block from (1, 1) to (n-nw, n-nw). The remaining
%       values corresponding to the covariances of the input 
%       disturbance model states are over-written at initialization.
%   R : (ny, ny) double
%       Output measurement noise covariance matrix.
%   nh : integer double
%       Number of hypotheses to model (each with a separate
%       Kalman filter).
%   n_min : integer double
%       Minimum life of cloned filters in number of sample
%       periods.
%   label : String (optional, default "MKF_SP_RODD")
%       Arbitrary name to identify observer instance.
%   d : integer double (optional, default 1)
%      Detection interval length in number of sample periods.
%   x0 : (n, 1) double (optional, default zeros)
%       Initial state estimates.
%   r0 : (1, 1) or (nh, 1) integer (optional)
%       Integer scalar or vector with values in the range 
%       {1, ..., nj} which indicate the prior system modes at time 
%       k = -1. If not provided, the default initialization based 
%       on the mode sequence that is generated will be used.
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

classdef MKFObserverSP_RODD < MKFObserverSP_DI
    properties (SetAccess = immutable)
        io struct
        nw double {mustBeInteger, mustBeNonnegative}
        n_shocks double {mustBeInteger, mustBeNonnegative}
    end
    properties
        sys_model struct
        %seq (:, 1) cell
        %i (1, 1) {mustBeInteger, mustBeNonnegative}
        %i_next (1, 1) {mustBeInteger, mustBeNonnegative}
        p_seq double
        Q0 double
        R double
        epsilon double  % TODO: make into cell like sigma_wp
        sigma_wp cell
    end
    methods
        function obj = MKFObserverSP_RODD(model,io,P0,epsilon, ...
                sigma_wp,Q0,R,nh,n_min,label,d,x0,r0)
            arguments
                model struct
                io struct
                P0 double
                epsilon double
                sigma_wp (1, :) cell
                Q0 (:, :) double
                R (:, :) double
                nh (1, 1) double {mustBeInteger}
                n_min (1, 1) double {mustBeInteger}
                label (1, 1) string = ""
                d (1, 1) double {mustBeInteger, ...
                    mustBeGreaterThanOrEqual(d,1)} = 1
                x0 = []
                r0 (:, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(r0, 1)} = 1
            end

            % Get model dimensions
            [n, nu, ny, ~, direct] = check_model(model);
            if direct
                assert(isequal(model.D, zeros(ny, nu)), ...
                    "ValueError: direct transmission not implemented")
            end
            assert(isequal(size(io.u_known), [nu, 1]))
            assert(isequal(size(io.y_meas), [ny, 1]))

            % Check size of initial process covariance matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Number of manipulatable inputs and unmeasured disturbances
            nu = sum(io.u_known);

            % Number of unmeasured disturbances
            nw = sum(~io.u_known);
            assert(nw > 0, "ValueError: io.u_known");

            % Number of switching input signals
            n_shocks = sum(cellfun(@(s) size(s, 2) > 1, sigma_wp));

            % Construct observer model without unmeasured disturbance
            % inputs
            Bu = model.B(:, io.u_known);
            if direct
                Du = model.D(:, io.u_known);
            else
                Du = zeros(ny,nu);
            end

            % Construct process noise covariance matrices for each 
            % possible input disturbance (returns a cell array)
            [Q, p_rk_g_rkm1] = construct_Q_model_SP(Q0, model.B, ...
                io.u_known, epsilon, sigma_wp);

            % Number of models (each with a hypothesis sequence)
            nj = numel(Q);

            % Transition probability matrix
            % Note that for RODD disturbances Pr(r(k)) is
            % assumed to be an independent random variable.
            T = repmat(p_rk_g_rkm1', nj, 1);

            % Initialize indicator sequences with 1s
            %seq = mat2cell(int16(ones(nh, nf)), int16(ones(1, nh)), nf);

            % System models are all the same - only Q switches
            model_obs = struct;
            model_obs.A = model.A;
            model_obs.B = Bu;
            model_obs.C = model.C;
            model_obs.D = Du;
            model_obs.Ts = model.Ts;
            models = repmat({model_obs}, 1, nj);
            for i = 1:nj
                models{i}.Q = Q{i};
                models{i}.R = R;
            end

            % Sequence pruning algorithm initialization
            % Assign all probability to first filter.
            p_seq_g_Yk_init = [1; zeros(nh-1, 1)];

            % Create MKF super-class observer instance
            obj = obj@MKFObserverSP_DI(models,P0,T,d,nh,n_min,label,x0, ...
                r0,p_seq_g_Yk_init,false);

            % Add additional variables used by this observer
            obj.sys_model = model;
            obj.io = io;
            obj.nw = nw;
            obj.n_shocks = n_shocks;
            obj.P0 = P0;
            obj.Q0 = Q0;
            obj.R = R;
            obj.epsilon = epsilon;
            obj.sigma_wp = sigma_wp;
            obj.type = "MKF_SP_RODD";

            % Initialize variables
            obj.reset()

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Call reset method of super class object
            reset@MKFObserverSP_DI(obj);

            % Set estimate covariances to high values for all the
            % filters except the first
            for i = 2:obj.nh
                obj.filters.Pkp1(:,:,i) = 1e10 * eye(obj.n);
            end

            % Reset sequence index
            %obj.i = int16(0);
            %obj.i_next = int16(1);

        end
    end
end
