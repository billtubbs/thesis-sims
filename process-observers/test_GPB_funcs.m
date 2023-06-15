% Test the following GPB functions:
%  - prob_transitions.m
%  - merge_estimates.m
%  - weighted_avg_estimates.m
%  - mode_transitions_all.m
%

clear all
rng(0)


%% Test prob_transitions.m

% Transition probability matrix (2 modes)
T = [0.95 0.05; 0.01 0.99];

% Transitions (x4)
rkm1 = [1 2 1 2]';
rk = [1 1 2 2]';
n_filt = size(rk, 1);

% Calculate transition probabilities
p_gamma_k_g_rkm1 = prob_transitions(rk, rkm1, T);
assert(isequal(p_gamma_k_g_rkm1, reshape(T, [], 1)));

% Simple way to calculate them
p_test = zeros(n_filt, 1);
for i = 1:n_filt
    p_test(i) = T(rkm1(i), rk(i));
end
assert(isequal(p_gamma_k_g_rkm1, p_test))

% Alternative way to calculate them
idx = sub2ind(size(T), rkm1, rk);
assert(isequal(idx, (1:4)'))
assert(isequal(p_gamma_k_g_rkm1, T(idx)))


%% Test merge_estimates.m

% Example 1
Xk = cat(3, 1, 1, 2, 2, 5);
Yk = cat(3, 2, 2, 4, 4, 10);
Pk = cat(3, 10, 10, 20, 20, 50);
p_seq_g_Yk = ones(5, 1) / 5;

% Simple tests with number of states, n = 1
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk);
assert(Xk_m == 2.2)
assert(Yk_m == 4.4)
assert(p_seq_g_Yk_m == 1)
assert(abs(Pk_m - 24.16) < 5e-15)

% Partial merge based on merge index
rk = [1 1 2 2 1]';
idx = [1 1 2 2 3]';
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);
assert(all(abs(Xk_m - cat(3, 1, 2, 5)) < 1e-15))
assert(all(abs(Pk_m - cat(3, 10, 20, 50)) < 1e-15))
assert(all(abs(Yk_m - cat(3, 2, 4, 10)) < 1e-15))
assert(all(abs(p_seq_g_Yk_m - [0.4 0.4 0.2]') < 1e-15))

% Same example but also returning merged modes rk_m
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m, rk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);
assert(all(abs(Xk_m - cat(3, 1, 2, 5)) < 1e-15))
assert(all(abs(Pk_m - cat(3, 10, 20, 50)) < 1e-15))
assert(all(abs(Yk_m - cat(3, 2, 4, 10)) < 1e-15))
assert(all(abs(p_seq_g_Yk_m - [0.4 0.4 0.2]') < 1e-15))
assert(isequal(rk_m, [1, 2, 1]'))

% Example 2
rk = [1 1 2 2 2]';
idx = [1 1 2 2 2]';
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);
assert(all(abs(Xk_m - cat(3, 1, 3)) < 1e-15))
assert(all(abs(Pk_m - cat(3, 10, 32)) < 1e-15))
assert(all(abs(Yk_m - cat(3, 2, 6)) < 1e-15))
assert(all(abs(p_seq_g_Yk_m - [0.4 0.6]') < 1e-15))

% 3
rk = [2 2 1 1 1]';
idx = [2 2 1 1 1]';
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m, rk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);
assert(all(abs(Xk_m - cat(3, 3, 1)) < 1e-15))
assert(all(abs(Pk_m - cat(3, 32, 10)) < 1e-15))
assert(all(abs(Yk_m - cat(3, 6, 2)) < 1e-15))
assert(all(abs(p_seq_g_Yk_m - [0.6 0.4]') < 1e-15))
assert(isequal(rk_m, [1 2]'))

% 4
p_seq_g_Yk = [2 2 1 1 0.4]' / 6.4;
rk = [1 1 2 2 3]';
idx = [1 1 2 2 3]';
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m, rk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);
assert(all(abs(Xk_m - cat(3, 1, 2, 5)) < 1e-15))
assert(all(abs(Pk_m - cat(3, 10, 20, 50)) < 1e-15))
assert(all(abs(Yk_m - cat(3, 2, 4, 10)) < 1e-15))
assert(all(abs(p_seq_g_Yk_m - [0.625 0.3125 0.0625]') < 1e-15))
assert(isequal(rk_m, [1, 2, 3]'))

% 5
p_seq_g_Yk = [0.2 0.2 0.2 0.2 0.2]';
rk = [1 2 2 3 3]';
idx = [1 2 2 3 3]';
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);
assert(all(abs(Xk_m - cat(3, 1, 1.5, 3.5)) < 1e-15))
assert(all(abs(Pk_m - cat(3, 10, 15.25, 37.25)) < 2e-15))
assert(all(abs(Yk_m - cat(3, 2, 3, 7)) < 1e-15))
assert(all(abs(p_seq_g_Yk_m - [0.2 0.4 0.4]') < 1e-15))

% 6
p_seq_g_Yk = [1 1 1 1 1]' / 5;
rk = [1 1 1 1 1]';
idx = [1 1 1 1 1]';
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);
assert(Xk_m == mean(Xk))
assert(abs(Pk_m - 24.16) < 1e-12)
assert(Yk_m == mean(Yk))
assert(p_seq_g_Yk_m == 1)

% With n = 2
% 1
x = [1 5]';
n = size(x,1);
Xk = cat(3, x, x+1, x+2, x+3);
A = reshape(1:4, 2, 2);
Pk = cat(3, A, 2.*A, 3.*A, 4.*A);
nh = size(Pk, 3);

rk = [1 1 2 2]';  % merging pattern
p_seq_g_Yk = [0.3 0.2 0.1 0.4]';
nj = 2;

% 1. Single state and output
Xk = cat(3, 1.1, 1.2, 1.4, 5);
Pk = cat(3, 100, 110, 120, 130);
Yk = cat(3, -2.2, -2.4, -2.8, -10);
p_seq_g_Yk = [0.2 0.2 0.5 0.1]';
rk = [1 1 2 2]';
idx = [1 1 2 2]';

% Function to test
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);

% Test outputs
Xk_m_test = cat(3, (0.2*1.1+0.2*1.2)/0.4, (0.5*1.4+0.1*5)/0.6);
assert(all(abs(Xk_m - Xk_m_test) < 1e-15))
assert(all(abs(Yk_m  - (-2) * Xk_m_test) < 2e-15))
Pk_m_test = cat(3, 0, 0);
Pk_m_test(:,:,1) = ...
      p_seq_g_Yk(1) .* (Pk(:,:,1) + (1.1 - Xk_m(:,:,1)).^2) / 0.4 ...
    + p_seq_g_Yk(2) .* (Pk(:,:,2) + (1.2 - Xk_m(:,:,1)).^2) / 0.4;
Pk_m_test(:,:,2) = ...
      p_seq_g_Yk(3) .* (Pk(:,:,3) + (1.4 - Xk_m(:,:,2)).^2) / 0.6 ...
    + p_seq_g_Yk(4) .* (Pk(:,:,4) + (5.0 - Xk_m(:,:,2)).^2) / 0.6;
assert(all(abs(Pk_m - Pk_m_test) < 2e-14))
assert(all(abs(p_seq_g_Yk_m - [0.4 0.6]') < 1e-15))

% 2. Multiple states
xk_est = [
    1.1 1.2 1.3 1.4; ...
    2.1 2.2 2.3 2.4; ...
    3.1 3.2 3.3 3.4 ...
];
n = size(xk_est, 1);
Xk = reshape(xk_est, 3, 1, 4);
Pk = cat(3, 10*eye(n), 20*eye(n), 30*eye(n), 40*eye(n));
p_seq_g_Yk = [0.3 0.2 0.1 0.4]';
rk = [1 1 2 2]';
idx = [1 1 2 2]';

% Function to test
[Xk_m, Pk_m, Yk_m, p_seq_g_Yk_m] = ...
    merge_estimates(Xk, Pk, Yk, p_seq_g_Yk, rk, idx);

% Compare with this code adapted from GPB2_estimation.m
Mix_u = reshape(p_seq_g_Yk, nj, 2);
for j = 1:nj

    % Normalize the mode probabilities
    Updx = xk_est(:, rk == j);
    Upd_u = Mix_u(:, j) / sum(Mix_u(:, j));
    UpdP = Pk(:,:, rk == j);

    % Mix the estimates
    Mix_x(:,j) = sum(Updx .* repmat(Upd_u', n, 1), 2);
    Mix_P(:,:,j) = zeros(n, n);
    for i = 1:nj
        summP = Upd_u(i) * (UpdP(:,:,i) + ...
            (Updx(:,i) - Mix_x(:,j)) * (Updx(:,i) - Mix_x(:,j))');
        Mix_P(:,:,j) =  Mix_P(:,:,j) + summP;
    end

    Storu(j) = sum(Mix_u(:,j));

end

% Tests
assert(isequal(Mix_x(:,1), ...
    sum(xk_est(:, 1:2) .* p_seq_g_Yk(1:2)', 2) / sum(p_seq_g_Yk(1:2))))
assert(isequal(Mix_x(:,2), ...
    sum(xk_est(:, 3:4) .* p_seq_g_Yk(3:4)', 2) / sum(p_seq_g_Yk(3:4))))
% Above is equivalent to
assert(isequal(Mix_x(:,1), ...
    sum(xk_est(:, rk == 1) .* p_seq_g_Yk(rk == 1)', 2) / ...
        sum(p_seq_g_Yk(rk == 1))))
assert(isequal(Mix_x(:,2), ...
    sum(xk_est(:, rk == 2) .* p_seq_g_Yk(rk == 2)', 2) / ...
        sum(p_seq_g_Yk(rk == 2))))

assert(all(abs(reshape(Xk_m, 3, 2) - Mix_x) < 1e-15, [1 2]));
assert(all(abs(Pk_m - Mix_P) < 1e-15, [1 2 3]));


%% Test weighted_avg_estimates.m
%
% TODO: Is this still used?  Or replaced by merge_estimates?

nj = 3;
n = 4;
ny = 2;
Xkf_est = randn(n, 1, nj);
Ykf_est = randn(ny, 1, nj);
Pkf_est = randn(n, n, nj);

p_seq_g_Yk = randn(nj, 1);
p_seq_g_Yk = p_seq_g_Yk ./ sum(p_seq_g_Yk);

[xk_est, yk_est, Pk] = weighted_avg_estimates(Xkf_est, Ykf_est, ...
    Pkf_est, p_seq_g_Yk);

assert(isequal(xk_est, Xkf_est(:,:,1) * p_seq_g_Yk(1) + ...
    Xkf_est(:,:,2) * p_seq_g_Yk(2) + Xkf_est(:,:,3) * p_seq_g_Yk(3)))
assert(isequal(yk_est, Ykf_est(:,:,1) * p_seq_g_Yk(1) + ...
    Ykf_est(:,:,2) * p_seq_g_Yk(2) + Ykf_est(:,:,3) * p_seq_g_Yk(3)))

Pk_m_test = zeros(n, n);
for i = 1:nj
    summP = p_seq_g_Yk(i) * (Pkf_est(:,:,i) + (Xkf_est(:,i) - xk_est) * (Xkf_est(:,i) - xk_est)');
    Pk_m_test = Pk_m_test + summP;
end
assert(isequal(Pk, Pk_m_test))


%% Test mode_transitions_all.m

[rkm1, rk] = mode_transitions_all(2);
assert(isequal(rkm1, [1 1 2 2]'))
assert(isequal(rk, [1 2 1 2]'))

[rkm1, rk] = mode_transitions_all(3);
assert(isequal(rkm1, [1 1 1 2 2 2 3 3 3]'))
assert(isequal(rk, [1 2 3 1 2 3 1 2 3]'))
