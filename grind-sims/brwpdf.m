function p = brwpdf(x, alpha1, alpha2, beta, tau, sd_e)
% Stationary probability density of the bounded random 
% walk (BRW). See Eq. 1 in Nicolau (2002)).
% Note: beta is the 'k' parameter in Nicolau's paper.
%
    p = sd_e^(-2) * exp( ...
        -2 * exp(beta) / sd_e^2 * (exp(-alpha1*(x - tau)) ...
        / alpha1 + exp(alpha2*(x - tau)) / alpha2) ...
    );
end