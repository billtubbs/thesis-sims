function [n, nu, ny] = check_dimensions(A, B, C, D)
% [n, nu, ny] = check_dimensions(A, B, C, D)
% Checks and returns dimensions of state-space system
% represented by matrices A, B, C, D.
%
    n = size(A, 1);
    assert(size(A, 2) == n, "ValueError: size(A)")
    nu = size(B, 2);
    assert(size(B, 1) == n, "ValueError: size(B)")
    ny = size(C, 1);
    assert(size(C, 2) == n, "ValueError: size(C)")
    assert(isequal(size(D), [ny nu]), "ValueError: size(D)")
end