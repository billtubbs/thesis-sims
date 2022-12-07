function p = prob_gamma_given_w(gamma_k, wk, epsilon, sigma_w)
% Probability of gamma(k) given w(k) calculated
% using Bayes rule.
%
% Example:
% >> epsilon = [0.99; 0.01];
% >> sigma_w = [0.01; 1];
% >> prob_gamma_given_w(0, 0.1, epsilon, sigma_w)
% 
% ans =
% 
%    1.9190e-18
% 

    p = prob_w_given_gamma(wk, gamma_k, sigma_w) ...
        * prob_gamma(gamma_k, epsilon) ./ prob_w(wk, epsilon, sigma_w);
return