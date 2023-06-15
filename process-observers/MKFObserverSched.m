% Multi-model Kalman Filter class definition
%
% Object class for simulating a multi-model Kalman filter for 
% state estimation with prescribed schedule to determine which 
% system is active at each time step. Although it is a multi-
% model observer because it has more than one possible system 
% model, it only uses one Kalman Filter (n_filt = 1).
%
% This scheduled KF is useful as a benchmark to compare other 
% multi-model observers with.
%

% TODO: Could this be inherited from MKFObserver?
%   or vice versa? SKF is an MKF without probability updates.

classdef MKFObserverSched < matlab.mixin.Copyable
    properties (SetAccess = immutable)
        Ts (1, 1) double {mustBeNonnegative}
        n (1, 1) double {mustBeInteger, mustBeNonnegative}
        nu (1, 1) double {mustBeInteger, mustBeNonnegative}
        ny (1, 1) double {mustBeInteger, mustBeNonnegative}
        nj (1, 1) double {mustBeInteger, mustBeNonnegative}
        d (1, 1) double {mustBeInteger, mustBeNonnegative}
        f (1, 1) double {mustBeInteger, mustBeNonnegative}
        n_filt (1, 1) double {mustBeInteger, mustBeNonnegative}
    end
    properties
        A cell
        B cell
        C cell
        D cell
        P0 double
        Q cell
        R cell
        K double
        Pkp1 double
        seq (1, :) double
        label (1, 1) string
        x0 (:, 1) double
        i (1, 2) {mustBeInteger, mustBeNonnegative}
        i_next (1, 2) {mustBeInteger, mustBeNonnegative}
        gamma_k double
        filter  KalmanFilter
        xkp1_est (:, 1) double
        ykp1_est (:, 1) double
        type (1, 1) string
    end
    methods
        function obj = MKFObserverSched(A,B,C,Ts,P0,Q,R,seq,label,x0)
        % obs = MKFObserverSched(A,B,C,Ts,P0,Q,R,seq,label,x0)
        %
        % Arguments:
        %	A, B, C : cell arrays
        %       discrete-time system matrices for each switching
        %       system modelled.
        %   Ts : double
        %       Sampling period.
        %   P0 : matrix, size (n, n)
        %       Initial value of covariance matrix of the state
        %       estimates.
        %   Q : cell array
        %       Process noise covariance matrices for each switching
        %       system.
        %   R : cell array
        %       Output measurement noise covariance matrices for each
        %       switching system.
        %   seq : row vector, size (1, f)
        %       Model indicator sequence. If there is only one random
        %       random input variable.
        %   d : detection interval length in number of sample periods.
        %   label : string (optional)
        %       Name.
        %   x0 : vector, size(n, 1), (optional)
        %       Intial state estimates.
        %

            % System dimensions
            [n, nu, ny] = check_dimensions(A{1}, B{1}, C{1});

            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.Ts = Ts;
            obj.P0 = P0;
            obj.Q = Q;
            obj.R = R;
            obj.seq = seq;
            obj.label = label;

            if nargin < 10
                x0 = zeros(n,1);
            end
            obj.x0 = x0;

            % Number of switching systems
            obj.nj = numel(A);

            % Check all other system matrix dimensions have same 
            % input/output dimensions and number of states.
            for j = 2:obj.nj
                [n_j, nu_j, ny_j] = check_dimensions(A{j}, B{j}, C{j});
                assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
                    "ValueError: size of A, B, and C")
            end

            % Initialize Kalman filter
            obj.filter = KalmanFilter(A{1},B{1},C{1},Ts,P0,Q{1},R{1}, ...
                "KF",x0);

            % System model indicator sequence
            obj.seq = seq;

            % Detection interval (not currently used)
            obj.d = 1;

            % Number of filters required (only 1)
            obj.n_filt = 1;

            % Fusion horizon length (full length of sequence)
            obj.f = size(obj.seq, 2);

            % Initialize all variables
            obj.reset()

            % Store other parameters
            obj.Ts = obj.filter.Ts;
            obj.n = obj.filter.n;
            obj.nu = obj.filter.nu;
            obj.ny = obj.filter.ny;
            obj.type = "SKF";

        end
        function reset(obj)
        % obj.reset()
        % Initialize all variables to initial values specified
        % when observer object was created.
        %

            % Switching variable at previous time instant
            obj.gamma_k = 0;

            % Initialize sequence index and counter
            % obj.i(1) is the sequence index (1 <= i(1) <= obj.f)
            % obj.i(2) is the counter for prob. updates (1 <= i(2) <= obj.d)
            obj.i = int16([0 0]);
            obj.i_next = int16([1 1]);

            % Initialize Kalman filter
            obj.filter.reset()

            % Initialize estimates
            obj.xkp1_est = obj.filter.xkp1_est;
            obj.ykp1_est = obj.filter.ykp1_est;

            % Initialize covariance matrix
            obj.Pkp1 = obj.filter.Pkp1;

        end
        function update(obj, yk, uk)
        % obj.update(yk, uk)
        % updates the scheduled Kalman filter and calculates the
        % estimates of the states and output at the next sample
        % time.
        %
        % Arguments:
        %   obs : struct containing the multi-model Kalman filter
        %       variables (see function mkf_filter).
        %   uk : vector (nu, 1) of system inputs at the current 
        %       sample time.
        %   yk : vector (ny, 1) of system output measurements
        %       at the current sample time.
        %

            assert(isequal(size(uk), [obj.nu 1]), "ValueError: size(uk)")
            assert(isequal(size(yk), [obj.ny 1]), "ValueError: size(yk)")

            % Increment sequence index and update counter
            % obj.i(1) is the sequence index (1 <= i(1) <= obj.f)
            % obj.i(2) is the counter for prob. updates (1 <= i(2) <= obj.d)
            % Whenever obj.i(2) exceeds obj.d (the spacing parameter), it is
            % reset to 1, and the sequence index obj.i(1) is incremented.
            obj.i = obj.i_next;
            obj.i_next = [mod(obj.i(1) - 1 + ...
                          idivide(obj.i(2), obj.d), obj.f) + 1, ...
                          mod(obj.i(2), obj.d) + 1];
            if obj.i(1) > obj.f
                error("Reached end of SKF sequence")
            end

            % Set system indicator based on provided schedule
            obj.gamma_k = obj.seq(obj.i(1));

            % Model index at current sample time
            ind = obj.gamma_k + 1;  % MATLAB indexing

            % Note: Currently, this update must happen every sample 
            % period so that S-functions do not have to memorize all
            % the model parameters each timestep.
            %if obj.i(2) == 1
            % Select filter system model based on current
            % model indicator value
            % TODO: only need to change these if ind changed?
            obj.filter.A = obj.A{ind};
            obj.filter.B = obj.B{ind};
            obj.filter.C = obj.C{ind};
            obj.filter.Q = obj.Q{ind};
            obj.filter.R = obj.R{ind};

            % Update observer estimates, gain and covariance matrix
            obj.filter.update(yk, uk);

            % Copy gain and covariance matrix
            obj.K = obj.filter.K;
            obj.Pkp1 = obj.filter.Pkp1;

            % Copy filter estimates
            obj.xkp1_est = obj.filter.xkp1_est;
            obj.ykp1_est = obj.filter.ykp1_est;

        end
    end
end
