function p = prob_rk(rk, p_rk)
% p = prob_rk(rk, p) returns the probability
% of the integer value or values of r(k) given the 
% probabilities in p which are defined as follows:
%   Pr(r(k)=1) = p(1)
%   Pr(r(k)=2) = p(2)
%   ...etc.
%
% For multiple independent discrete variables, define 
% probabilities  in the columns of epsilon and sequence 
% values in the rows of r(k).
%
% Example 1 - single binary variable (two states).
% >> prob_rk(2, [0.99; 0.01])
% 
% ans =
% 
%     0.0100
%
% Example 2 - a sequence of two discrete variables with three states
% >> prob_rk([2 1 1 1; 1 3 1 1], [0.98 0.95; 0.01 0.03; 0.01 0.02])
% 
% ans =
%
%     0.0100    0.9800    0.9800    0.9800
%     0.9500    0.0200    0.9500    0.9500
% 

    n = size(rk, 1);
    assert(n == size(p_rk, 2))
    p = zeros(size(rk));
    for i = 1:n
        p(i, :) = p_rk(rk(i, :), i);
    end

end