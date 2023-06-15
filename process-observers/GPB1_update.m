function [xk_est,yk_est,Pk,p_seq_g_Yk] = GPB1_update(models,T, ...
    Xkp1f_est,Pkp1f,yk,p_seq_g_Yk)
% [xk_est,yk_est,Pk,p_seq_g_Yk] = GPB1_update(models,T,Xkp1f_est, ...
%     Pkp1f,yk,p_seq_g_Yk)
% Update equations for simulating the first-order 
% generalised pseudo-Bayes (GPB1) multi-model Kalman 
% filter for state estimation of Markov jump linear
% systems.
%
% Based on code from the following article:
% -  Bayesian State Estimations for Markovian Jump Systems, 
%    by Shunyi Zhao, Choon Ki Ahn, Peng Shi, Yuriy S. Shmaliy, 
%    and Fei Liu, 2019.
%

    % TODO: better to move this code into the main object class
    % as a method?  (filters.Kf and filters.Sk are not stored).

    % All models must have the same dimensions (this is not checked)
    nj = numel(models);  % number of models (= no. of system modes)
    n = size(models{1}.A, 1);
    ny = size(models{1}.C, 1);

    updx = nan(n, nj);
    updy = nan(ny, nj);
    updP = nan(n, n, nj);
    p_yk_g_seq_Ykm1 = nan(nj, 1);

    for j = 1:nj

        xkp1_est = Xkp1f_est(:,:,j);
        Pkp1 = Pkp1f(:,:,j);

        % Select model
        m = models{j};

        % Kalman filter correction step
        [updx(:,j),updP(:,:,j),updy(:,j),Kf,Sk] = ...
            kalman_update_f(m.C,m.R,xkp1_est,Pkp1,yk);

        % Output estimate
        ykp1_est = m.C * xkp1_est;

        % Update the conditional likelihood given the data
        p_yk_g_seq_Ykm1(j) = mvnpdf(yk, ykp1_est, Sk);

    end

    % Compute likelihoods of models given data
    p_seq_g_Ykm1 = sum(T .* p_seq_g_Yk, 1)';
    cond_pds = p_seq_g_Ykm1 .* p_yk_g_seq_Ykm1;
    p_seq_g_Yk = cond_pds / sum(cond_pds);

    % Compute updated merged estimates
    xk_est = updx * p_seq_g_Yk;
    % Above is equivalent to:
    %   sum(updx .* repmat(p_seq_g_Yk',n,1), 2)
    yk_est = updy * p_seq_g_Yk;
    % Above is equivalent to:
    %   sum(updy .* repmat(p_seq_g_Yk',ny,1), 2)
    Pk = zeros(n,n);
    for i = 1:nj
        Pk = Pk + p_seq_g_Yk(i) * (updP(:,:,i) + ...
            (updx(:,i) - xk_est) * (updx(:,i) - xk_est)');
    end

end
