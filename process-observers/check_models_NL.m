function [nj, n, nu, ny, Ts] = check_models_NL(models)
% [nj, n, nu, ny, Ts] = check_models_NL(models)
% Checks all system models have same input/output 
% dimensions and number of states and returns the
% dimensions.
%
% Arguments:
%   models : (1, nj) cell array
%       Cell array of structs containing the functions and 
%       parameters of the non-linear models of the switching
%       system. Each model struct must include:
%           model.f - the transition function,
%           model.h - the measurement function,
%           model.dfdx - partial derivative of f,
%           model.dhdx - partial derivative of f,
%           model.n - integer for the number of model 
%               states,
%           model.nu - integer for the number of inputs,
%           model.ny - integer for the number of outputs,
%           model.Ts - the sampling period,
%           model.params (optional) - cell array of fixed 
%               parameters that must be passed to f, h, 
%               dfdx, and dhdx.
%

    % Number of models
    nj = numel(models);

    % Get dimensions of first model
    [n, nu, ny, Ts] = check_model_NL(models{1});

    % Check other models have same dimensions and sampling 
    % time
    for j = 2:nj
        [n_j, nu_j, ny_j, Ts_j] = check_model_NL(models{j});
        assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
            "ValueError: mismatch in model dimensions")
        assert(Ts_j == Ts, "ValueError: Ts")
    end

end