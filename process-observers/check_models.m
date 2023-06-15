function [nj, n, nu, ny, Ts] = check_models(models)
% [nj, n, nu, ny, Ts] = check_models(models)
% Checks all system models have same input/output 
% dimensions and number of states and returns the
% dimensions.
%
% Arguments:
%   models : (1, nj) cell array of structs
%       Each struct contains the parameters of a linear
%       dynamical system model. These include: A, B, 
%       and C for the system matrices, Q and R for the
%       state error covariance and output measurement 
%       noise covariance, and Ts for the sample period.
%

    % Number of models
    nj = numel(models);

    % Get dimensions of first model and see if it has direct
    % transmission (D matrix)
    [n, nu, ny, Ts, direct] = check_model(models{1});

    % Check other models have same dimensions and sampling 
    % time
    for j = 2:nj
        [n_j, nu_j, ny_j, Ts_j, d_j] = check_model(models{j});
        assert(isequal([n_j, nu_j, ny_j], [n, nu, ny]), ...
            "ValueError: size of A, B, and C")
        assert(direct == d_j, "ValueError: D")
        assert(Ts_j == Ts, "ValueError: Ts")
    end

end