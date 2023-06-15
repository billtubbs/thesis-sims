function [n, nu, ny, Ts] = check_model_NL(model)
% [n, nu, ny, Ts] = check_model_NL(model)
% Get dimensions from model struct, check they are 
% consistent and determine if there is a D matrix 
% (direct transmission), and return the sampling
% period.
%
% Arguments:
%   model : struct
%       Struct containing functions and other parameters
%       of a non-linear model of the system. This must 
%       include:
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
% Returns:
%   n : integer
%       Number of states
%   nu : integer
%       Number of inputs
%   ny : integer
%       Number of outputs
%   Ts : double
%       Sampling period
%

    % Model dimensions
    n = model.n;
    nu = model.nu;
    ny = model.ny;

    % Sampling period
    Ts = model.Ts;

    % Check struct has additional required attributes
    assert(isa(model.f, 'function_handle'), ...
        "ValueError: model.f")
    assert(isa(model.h, 'function_handle'), ...
        "ValueError: model.h")

end