function [xk_est, yk_est, Pk] = weighted_avg_estimates(Xkf_est, ...
    Ykf_est, Pkf_est, p_seq_g_Yk)
% Compute multi-model observer state and output estimates
% and estimated state error covariance using the weighted-
% averages based on the conditional probabilities.
%
    weights = reshape(p_seq_g_Yk, 1, 1, []);
    xk_est = sum(weights .* Xkf_est, 3);
    yk_est = sum(weights .* Ykf_est, 3);
    Xkf_devs = xk_est - Xkf_est;
    Pk = sum(weights .* (Pkf_est + ...
        pagemtimes(Xkf_devs, pagetranspose(Xkf_devs))), 3);
end