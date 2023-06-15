function model = model_add_outdist(model, y_meas, y_scale)
% model = model_add_outdist(model, y_meas)
% Augment the discrete-time state-space model with
% input disturbances on each manipulatable input.
%

    % Get model dimensions
    [n, nu, ny, Ts] = check_model(model);

    if nargin < 3
        % Default: no scaling of outputs
        y_scale = ones(1, ny);
    end

    % Augmented model with integrators added to each
    % measured output
    Aaug = [model.A zeros(n,ny);
            zeros(ny,n) eye(ny)];
    Baug = [model.B; 
            zeros(ny,nu)];
    Caug = [model.C y_scale.*eye(ny)];
    Daug = [zeros(ny,nu)];
    model = ss(Aaug,Baug,Caug,Daug,Ts);

end