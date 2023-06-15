function [f,h,params] = make_trans_funcs_ssd(A,B,C,D)
% [f,h,params] = make_trans_funcs_ssd(A,B,C,D)
% Make state transition and measurement functions
% for a linear system from the state-space matrices
% (f, h, params are used by run_simulation_obs)
    params.A = A;
    params.B = B;
    params.C = C;
    params.D = D;
    f = @(xk, uk, dt, params) params.A * xk + params.B * uk;
    h = @(xk, uk, dt, params) params.C * xk + params.D * uk;
end

% TODO: These functions might need time t as argument.
% TODO: Make a unit test script