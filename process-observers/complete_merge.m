function [xk_est,Pk,yk_est] = complete_merge(Xkf_est,Pkf,Ykf_est,p_seq_g_Yk)
    weights = reshape(p_seq_g_Yk, 1, 1, []);
    xk_est = sum(weights .* Xkf_est, 3);
    yk_est = sum(weights .* Ykf_est, 3);
    Xkf_devs = xk_est - Xkf_est;
    Pk = sum(weights .* (Pkf + ...
        pagemtimes(Xkf_devs, pagetranspose(Xkf_devs))), 3);
end