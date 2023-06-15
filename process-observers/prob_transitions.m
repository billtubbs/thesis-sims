function p = prob_transitions(gamma_k, gamma_km1, T)
% p = prob_transitions(gamma_k, gamma_km1, T) returns the 
% probabilities of the mode transitions from gamma_km1 to
% gamma_k given the transition probabilitiy matrix T.
% defined as:
%   Pr(gamma(k)==j | gamma(k-1)==i) = T(i, j)
%
% Example 1 - transitions of a binary variable.
% >> T = [0.95 0.05; 0.01 0.99];
% >> gamma_km1 = [1 2 1 2]';
% >> gamma_k = [1 1 2 2]';
% >> prob_transitions(gamma_k, gamma_km1, T)
% 
% ans =
% 
%     0.9500
%     0.0100
%     0.0500
%     0.9500
%
    n = size(gamma_k, 1);
    assert(size(gamma_km1, 1) == n)
    p = zeros(n, 1);
    for i = 1:n
        p(i) = T(gamma_km1(i), gamma_k(i));
    end

end