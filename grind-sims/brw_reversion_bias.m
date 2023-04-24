function a = brw_reversion_bias(x, alpha1, alpha2, beta, tau)
% This is the function 'a(x)' from Nicolau (2002) used in the
% difference equation of the bounded random walk (BRW) 
% (see Eq. 1 in the paper).
% Note: beta is the 'k' parameter in Nicolau's paper.
%
    a = exp(beta) * (exp(-alpha1 * (x - tau)) - exp(alpha2 * (x - tau)));
end
