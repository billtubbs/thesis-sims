function [x, a] = sample_random_shocks(sz, p, sigma1, sigma2)
% [x, a] = sample_random_shocks(s, p, sigma1, sigma2)
%
% Generate an array of random shock values for one or more
% random variables x(k,i) as defined in section 2.1 of MacGregor
% et al. (1984) and Robertson et al. (1995).
%
% Arguments:
%   sz : array or scalar
%       the size of the output array [nT na] where nT is the
%       number of samples and na is the number of indepedent
%       random variables (na >= 1). If s is a scalar then it
%       is assumed that na = 1 and sz = nT.
%   p : scalar (if na == 1) or vector (na > 1)
%       p(i) is the probability of a shock occuring for x(k,i) 
%       (a(k,i) = 1) (0 < p < 1).
%   sigma1 : (optional) scalar (if na == 1) or vector (na > 1)
%       sigma1(i) are the standard deviations of the values of 
%       x(k,i) when a shock occurs (with probability p(i)). 
%       If not provided, sigma1(i) = 1.
%   sigma2 : (optional) scalar (if na == 1) or vector (na > 1)
%       sigma2(i) are the standard deviations of the values of 
%       x(k,i) when no shock occurs (with probability 1 - 
%       p(i)). If not provided, sigma2(i) = 0.
%
% Returns:
%   x : array
%       random disturbances [x(1); x(2); ... x(nT)].
%   a : array (optional)
%       binary values [a(1); a(2); ... a(nT)] indicating
%       when shocks occurred (0 = no shock, 1 = shock).
%
% Examples:
% >> sample_random_shocks(5, 0.2)
% 
% ans =
% 
%          0
%          0
%    -1.3499
%          0
%          0
% 
% >> [x, a] = sample_random_shocks([5 2], [0.8; 0.2], [10; 1])
% 
% x =
% 
%          0         0
%          0         0
%   -11.4707         0
%          0         0
%   -10.6887         0
% 
% 
% a =
% 
%   5×2 logical array
% 
%    0   0
%    0   0
%    1   0
%    0   0
%    1   0
% 
% References:
%  -  MacGregor, J. F., Harris, T. J., & Wright, J. D. (1984). 
%     Duality Between the Control of Processes Subject to Randomly 
%     Occurring Deterministic Disturbances and ARIMA Stochastic  
%     Disturbances. Technometrics, 26(4), 389-397.  
%     https://doi.org/10.1080/00401706.1984.10487992
%  -  Robertson, D. G., Kesavan, P., & Lee, J. H. (1995). 
%     Detection and estimation of randomly occurring 
%     deterministic disturbances. Proceedings of 1995 American
%     Control Conference - ACC'95, 6, 4453-4457. 
%     https://doi.org/10.1109/ACC.1995.532779
%

    if isscalar(sz)
        sz = [sz 1];
    end
    if nargin < 3
        sigma1 = ones(sz(2), 1);
    end
    if nargin < 4
        sigma2 = zeros(sz(2), 1);
    end
    if isscalar(sigma1)
        sigma1 = sigma1*ones(sz(2), 1);
    end
    if isscalar(sigma2)
        sigma2 = sigma2*ones(sz(2), 1);
    end
    a = rand(sz) < p';
    x = repmat(sigma2', sz(1), 1) .* randn(sz);
    for i=1:size(x, 2)
        n = sum(a(:, i));
        x(a(:,i), i) = sigma1(i) * randn(n, 1);
    end
end
