function p = prob_gamma(gamma, epsilon)
% p = prob_gamma(gamma, epsilon) returns the probability
% of the integer value or values of gamma given the probabilities
% in epsilon which are defined as follows:
%   Pr(gamma=0) = epsilon(1)
%   Pr(gamma=1) = epsilon(2)
%   ...etc.
% Examples 1 gamma is a scalar and epsilon is a vector.
% >> prob_gamma(1, [0.99; 0.01])
% 
% ans =
% 
%     0.0100
%
% Example 2 - gamma and epsilon are matrices
% >> prob_gamma([1 0 0; 0 1 0], [0.99 0.95; 0.01 0.05])
% 
% ans =
% 
%     0.0100    0.9900    0.9900
%     0.9500    0.0500    0.9500
% 

    assert(size(gamma, 1) == size(epsilon, 2))
    p = zeros(size(gamma));
    for i = 1:size(gamma, 1)
        p(i, :) = epsilon(gamma(i, :) + 1 + (i-1)*2);
    end

end