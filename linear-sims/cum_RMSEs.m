
function crmse = cum_RMSEs(E)
% cm = cumrmse(errors)
% Computes the root-mean-squared-error at times 1 to k
% for all k in 1:size(errors, 1).  E should be
% a column vector or matrix of errors computed at each
% timestep.
    n = size(E, 1);
    crmse = sqrt(cumsum(E.^2, 1) ./ (1:n)');
end