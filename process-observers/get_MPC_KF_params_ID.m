function [Q,R,Gpred,N] = get_MPC_KF_params_ID(mpcobj)
% [Q, R, Gpred] = get_MPC_KF_params_ID(mpcobj)
% Determines the parameters of the Kalman filter of the
% MPC object in the case of an input disturbance
% model. Also returns the discrete-time prediction 
% model used by the estimator.
%

    % Get scaling factors on MVs, DVs, and CVs
    MV_scale_factors = extractfield(mpcobj.ManipulatedVariables,'ScaleFactor');
    %DV_scale_factors = extractfield(mpcobj.DisturbanceVariables,'ScaleFactor');
    OV_scale_factors = extractfield(mpcobj.OutputVariables,'ScaleFactor');

    % Get intput disturbance model
    Gid = getindist(mpcobj);

    % Get measurement noise model
    if ~isempty(mpcobj.Model.Noise)
        Gmn = mpcobj.Model.Noise;
    else
        % Default measurement noise model
        ny = size(mpcobj,'mo');
        Gmn = ss(eye(ny).*OV_scale_factors,'Ts',mpcobj.Ts);
    end

    % Make sure there are no output disturbances
    assert(isempty(getoutdist(mpcobj)))

    % Get plant model
    if isct(mpcobj.Model.Plant)
        % Convert to discrete time model
        Gd = c2d(mpcobj.Model.Plant, mpcobj.Ts);
    else
        Gd = mpcobj.Model.Plant;
    end

    % Get dimensions of system model
    [n, ~, ny, Ts] = check_model(Gd);
    MVs = Gd.InputGroup.Manipulated;
    nu = numel(MVs);
    UDs = Gd.InputGroup.Unmeasured;
    nw = numel(UDs);

    % Augment the plant model with the disturbance model
    A = [Gd.A Gd.B(:,UDs)*Gid.C; zeros(size(Gid.A,1),size(Gd.A,2)) Gid.A];
    Bu = [Gd.B(:,MVs); zeros(ny,nu)];
    Cm = [Gd.C zeros(ny,nw)];
    D = Gd.D(:,MVs);
    % This should produce the same system as
    % Gpred = series([1 0; 0 Gid], Gd)

    % Prediction model 
    Gpred = ss(A,Bu,Cm,D,Ts);

    % Noise input equations
    %
    %   w[k] = [ Gd.B*Gid.D  0 ] * [ wn_id ] = B_est * white noise
    %          [ Gid.B       0 ]   [ wn_mn ]
    %
    %   v[k] = [ Gd.D Gmn.D ] * [ wn_id ] = D_est * white noise
    %                           [ wn_mn ]
    %

    % Noise input matrix - see equation for w(k)
    B_est = padarray(blkdiag(Gd.B(:,MVs).*MV_scale_factors, Gid.B), [0 size(Gmn.B,2)], 'post');
    %B_est = padarray([Gd.B(:,MVs).*MV_scale_factors; Gid.B], [0 size(Gmn.B,2)], 'post');

    % Noise direct transmission matrix - see equation for v(k)
    D_est = [Gd.D(:,MVs).*MV_scale_factors Gid.D Gmn.D];
    %D_est = [Gd.D(:,MVs).*MV_scale_factors Gmn.D];

    % Kalman filter noise covariance parameters
    Q = B_est * B_est';
    R = D_est * D_est';
    N = B_est * D_est';
    G = eye(n+nw);
    H = zeros(ny,n+nw);
    [~, L, ~, M] = kalman(ss(A,[Bu G],Cm,[D H],Ts),Q,R,N);

    % Check gains and model are identical to those of the MPC
    [L1,M1,A1,Cm1,Bu1] = getEstimator(mpcobj);
    assert(max(abs(L - L1), [], [1 2]) < 1e-10)
    assert(max(abs(M - M1), [], [1 2]) < 1e-10)
    assert(max(abs(A - A1), [], [1 2]) < 1e-10)
    assert(max(abs(Cm - Cm1), [], [1 2]) < 1e-10)
    assert(max(abs(Bu1 - Bu), [], [1 2]) < 1e-10)

end