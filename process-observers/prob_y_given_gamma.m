function p = prob_y_given_gamma(yk, gamma_k, uk, A, B, epsilon, sigma_w)
    xk = 
    p = prob_w_given_gamma(wk, gamma_k, sigma_w) ...
        * prob_gamma(gamma_k, epsilon) ./ prob_w(wk, sigma_w);
return