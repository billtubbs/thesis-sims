function model = model_add_indist(model, u_known)
% model = model_add_indist(model, u_known)
% Augment the discrete-time state-space model with
% input disturbances on each manipulatable input.
%

    % Get model dimensions
    [n, nu, ny, Ts] = check_model(model);

    % Augmented model with integrators added to each
    % known input
    Aaug = [model.A model.B(:,u_known); 
            zeros(nu,n) eye(nu)];
    Baug = [model.B(:,u_known);
            zeros(nu,nu)];
    Caug = [model.C zeros(ny,nu)];
    Daug = [zeros(ny,nu)];
    model = ss(Aaug,Baug,Caug,Daug,Ts);

end