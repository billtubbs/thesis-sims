function [n, nu, ny, Ts, direct] = check_model(model)
% [n, nu, ny, direct, Ts] = check_model(model)
% Get dimensions from model struct, check they are 
% consistent and determine if there is a D matrix 
% (direct transmission), and return the sampling
% period.
%
% Arguments:
%   model : struct
%       Struct containing the parameters of a linear
%       model of the system dynamics. These include: 
%       A, B, C, (and optionally, D), which are the 
%       system matrices, and the sampling period, Ts.
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
%   direct : logical
%       True if model has a D matrix.
%

    if isfield(model, "D") || isprop(model, "D")
        [n, nu, ny] = check_dimensions(model.A, model.B, model.C, model.D);
        direct = true;
    else
        [n, nu, ny] = check_dimensions(model.A, model.B, model.C);
        direct = false;
    end
    Ts = model.Ts;

end