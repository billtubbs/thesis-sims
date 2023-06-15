% Multi-model Kalman Filter class definition
%
% obs = MKFObserverSF_RODD95(model,io,P0,epsilon, ...
%     sigma_wp,Q0,R,f,m,d,label,x0,r0)
%
% Object class for simulating a multi-model observer for 
% state estimation in the presence of randomly-occurring 
% deterministic disturbances (RODDs) as described in 
% Robertson et al. (1995). This version is slightly 
% different to that described in Robertson et al. (1998).
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
%   f : integer double
%       Fusion horizon length in number of samples.
%   m : integer double
%       Maximum number of disturbances over fusion horizon.
%   d : integer double
%       Detection interval length in number of sample periods.
%   label : String (optional, default "MKF_SF_RODD")
%       Arbitrary name to identify observer instance.
%   r0 : (1, 1) or (nh, 1) integer (optional)
%       Integer scalar or vector with values in the range 
%       {1, ..., nj} which indicate the prior system modes at time 
%       k = -1. If not provided, the default initialization based 
%       on the mode sequence that is generated will be used.
%
% References:
%  -  Robertson, D. G., Kesavan, P., & Lee, J. H. (1995). 
%     Detection and estimation of randomly occurring 
%     deterministic disturbances. Proceedings of 1995 American
%     Control Conference - ACC'95, 6, 4453-4457. 
%     https://doi.org/10.1109/ACC.1995.532779
%  -  Robertson, D. G., & Lee, J. H. (1998). A method for the
%     estimation of infrequent abrupt changes in nonlinear 
%     systems. Automatica, 34(2), 261-270.
%     https://doi.org/10.1016/S0005-1098(97)00192-1%
%

classdef MKFObserverSF_RODD95 < MKFObserverSF
    properties (SetAccess = immutable)
        io struct
        nw double {mustBeInteger, mustBeNonnegative}
        n_shocks double {mustBeInteger, mustBeNonnegative}
        f double {mustBeInteger, mustBeNonnegative}
        m double {mustBeInteger, mustBeNonnegative}
        d (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        sys_model struct
        alpha double
        beta (1, 1) double
        p_seq double
        p_rk double  % TODO: Has this been replaced with p_rk_g_rkm1?
        Q0 {mustBeNumeric}
        R double
        epsilon double
        sigma_wp cell
    end
    methods
        function obj = MKFObserverSF_RODD95(model,io,P0,epsilon, ...
                sigma_wp,Q0,R,f,m,d,label,x0,r0)
            arguments
                model struct
                io struct
                P0 double
                epsilon double
                sigma_wp (1, :) cell
                Q0 (:, :) double
                R (:, :) double
                f (1, 1) double {mustBeInteger, mustBeGreaterThan(f, 0)}
                m (1, 1) double {mustBeInteger, mustBeGreaterThan(m, 0)}
                d (1, 1) double {mustBeInteger, mustBeGreaterThan(d, 0)}
                label (1, 1) string = ""
                x0 = []
                r0 (:, 1) int16 {mustBeGreaterThan(r0, 0)} = []
            end

            % Get model dimensions
            [n, nu, ny, ~, direct] = check_model(model);
            if direct
                assert(isequal(model.D, zeros(ny, nu)), ...
                    "ValueError: direct transmission not implemented")
            end
            assert(isequal(size(io.u_known), [nu, 1]), ...
                "ValueError: size(io.u_known)")
            assert(isequal(size(io.y_meas), [ny, 1]), ...
                "ValueError: size(io.y_meas)")

            % Check size of process covariance default matrix
            assert(isequal(size(Q0), [n n]), "ValueError: size(Q0)")

            % Number of manipulatable inputs and unmeasured disturbances
            nu = sum(io.u_known);

            % Number of unmeasured disturbances
            nw = sum(~io.u_known);
            assert(nw > 0, "ValueError: io.u_known");

            % Number of switching input signals
            n_shocks = sum(cellfun(@(s) size(s, 2) > 1, sigma_wp));
            assert(n_shocks == nw, ...
                "NotImplementedError: only switching disturbances.")

            % Construct observer model without unmeasured disturbance
            % inputs
            Bu = model.B(:, io.u_known);
            Bw = model.B(:, ~io.u_known);
            if direct
                Du = model.D(:, io.u_known);
            else
                Du = zeros(ny,nu);
            end

            % Probability of at least one shock in a detection interval
            % (Detection interval is d sample periods in length).
            if d == 1
                alpha = epsilon;
            else
                alpha = 1 - (1 - epsilon).^d;
            end

            % Convert fusion horizon to number of detection intervals
            assert(rem(f, d) == 0, "ValueError: Fusion horizon and " ...
                + "detection interval not compatible")
            nf = f / d;

            % Check for invalid parameter combinations
            assert(m <= nf, "ValueError: m <= f/d")

            % Construct process noise covariance matrices and switching
            % sequences over the fusion horizon, and the prior 
            % probabilities of each sequence.
            %TODO: Adapt this to allow non-switching random disturbances
            idx_switch = 1:nw;
            sigma_wp_RODDs = cell2mat(sigma_wp(idx_switch)');
            [Q, p_rk, S, p_seq] = construct_Q_model_SF95(Q0, Bw, alpha, ...
                sigma_wp_RODDs, nf, m, nw);

            % Number of models (each with a different hypothesis sequence)
            nj = numel(Q);

            % Number of hypotheses (filters) to be modelled
            nh = size(S, 1);

            % Expand sequences by inserting sample times with 'no-shock'
            % between times when shocks occur.
            seq = cell(nh, 1);
            for i = 1:nh
                seq{i} = int16(ones(size(S{i}, 1), f));
                % Add shock indications at start of each detection
                % interval
                seq{i}(:, 1:d:f) = S{i};
                % Alternatively, at end of each detection interval
                %seq{i}(:, d:d:f) = S{i};
            end
            assert(size(seq{1}, 2) == f)

            % Transition probability matrix
            % Note that for RODD disturbances Pr(gamma(k)) is
            % assumed to be an independent random variable.
            T = repmat(p_rk', nj, 1);

            % Tolerance parameter (total probability of defined sequences)
            beta = sum(p_seq);

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

            % Create MKF super-class observer instance
            obj = obj@MKFObserverSF(models,P0,seq,T,label,x0,r0,[],false);

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
            obj.f = f;
            obj.m = m;
            obj.d = d;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.p_seq = p_seq;
            obj.p_rk = p_rk;
            obj.type = "MKF_SF_RODD95";

            % Initialize variables
            obj.reset()

        end
    end
end
