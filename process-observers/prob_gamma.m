function p = prob_gamma(gamma, epsilon)
% p = prob_gamma(gamma, epsilon) returns the probability
% of the integer value or values of gamma given the 
% probabilities in epsilon which are defined as follows:
%   Pr(gamma=0) = epsilon(1)
%   Pr(gamma=1) = epsilon(2)
%   ...etc.
%
% For multiple independent discrete variables, define 
% probabilities  in the columns of epsilon and sequence 
% values in the rows of gamma.
%
% Example 1 - single binary variable (two states).
% >> prob_gamma(1, [0.99; 0.01])
% 
% ans =
% 
%     0.0100
%
% Example 2 - a sequence of two discrete variables with three states
% >> prob_gamma([1 0 0 0; 0 2 0 0], [0.98 0.95; 0.01 0.03; 0.01 0.02])
% 
% ans =
%
%     0.0100    0.9800    0.9800    0.9800
%     0.9500    0.0200    0.9500    0.9500
% 

    n = size(gamma, 1);
    assert(n == size(epsilon, 2))
    p = zeros(size(gamma));
    for i = 1:n
        p(i, :) = epsilon(gamma(i, :) + 1, i);
    end

end