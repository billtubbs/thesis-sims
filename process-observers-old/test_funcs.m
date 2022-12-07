% Test the following functions:
%
%  - prob_gamma.m
%  - prob_w_given_gamma.m
%  - prob_w.m
%  - prob_gamma_given_w.m
%  - n_filters.m
%  - seq_prob.m
%  - combinations.m 
%  - combinations_lte.m
%  - prob_combs.m
%

clear all

epsilon = [0.99; 0.01];
sigma1 = 0.01;
b = 100;
sigma_w = [sigma1; sigma1*b];


%% Test probability functions

% Test prob_gamma.m
assert(prob_gamma(0, epsilon) == 0.99)
assert(prob_gamma(1, epsilon) == 0.01)

% Test prob_w_given_gamma.m
assert(round(prob_w_given_gamma(0, 0, sigma_w), 4) == 39.8942);
assert(round(prob_w_given_gamma(0, 1, sigma_w), 4) == 0.3989);
assert(round(prob_w_given_gamma(1, 0, sigma_w), 4) == 0.0);

% Test prob_w.m
assert(round(prob_w(0, epsilon, sigma_w), 4) == 39.4993)
assert(all(round(prob_w([0.1; 1; -10], epsilon, sigma_w), 4) == [0.0040; 0.0024; 0.0000]))
assert(prob_w(0.01, epsilon, sigma_w) == prob_w_given_gamma(0.01, 0, sigma_w) * epsilon(1) ...
    + prob_w_given_gamma(0.01, 1, sigma_w) * epsilon(2))

% Test prob_gamma_given_w.m
wk = 0; gamma_k = 0;
assert(round(prob_gamma_given_w(wk, gamma_k, epsilon, sigma_w), 4) == 0.9999)
wk = 1; gamma_k = 1;
assert(round(prob_gamma_given_w(wk, gamma_k, epsilon, sigma_w), 4) == 1)

% Tests with vector inputs
w = linspace(-3, 3, 101);
p1 = normpdf(w,0,sigma_w(1));
p2 = normpdf(w,0,sigma_w(2));
assert(all(prob_w_given_gamma(w, 0, sigma_w) == p1));
assert(all(prob_w_given_gamma(w, 1, sigma_w) == p2));
assert(all(prob_w(w, epsilon, sigma_w) == (epsilon(2)*p2 + epsilon(1)*p1)))


%% Test n_filters.m and seq_prob.m

% Calculations with example settings from Robertson paper

nw = 2;  % dimension of w(k)
epsilon = 0.005;   % probability of shocks at each sample time
% (assume probability is same for all elements of w(k))
%assert(numel(epsilon) == nw);

d = 40;  % length of periods
p = (1 - epsilon).^d;  % probability of shock in each period
assert(round(p, 2) == 0.82)

m = 3;  % maximum number of shocks
n = 3;  % fusion horizon

assert(n_filters(m, n, nw) == 42)
assert(n_filters(m, n) == 8)

% Calculation method in Robertson et al.
p_seq = 0;
for k=0:m
    p_seq = p_seq + nchoosek(n*nw, k) * p^(n*nw - k) * (1 - p)^k;
end

assert(round(p_seq, 2) == 0.99)
assert(seq_prob(m, n*nw, p) == p_seq)


% Example 2 - single disturbance
nw = 1;
epsilon = 0.01;
d = 5;
m = 3;
n = 3;
n_filt = n_filters(m, n);
assert(n_filt == 8)

p = (1 - epsilon).^d;
p_seq = seq_prob(m, n*nw, p);
assert(p_seq == 1)


% Example 3 - single disturbance use for report simulation
nw = 1;
epsilon = 0.0100;
d = 5;
m = 2;
n = 6;
n_filt = n_filters(m, n);
assert(n_filt == 22)
p = (1 - epsilon).^d;

% Calculation method in Robertson et al.
p_seq = 0;
for k=0:m
    p_seq = p_seq + nchoosek(n*nw, k) * p^(n*nw - k) * (1 - p)^k;
end

assert(round(p_seq, 4) == 0.9979)
assert(seq_prob(m, n*nw, p) == p_seq)


%% Test combinations.m, combinations_lte.m, combinations_gte_prob
% and prob_Gammas.m

epsilon = [0.99; 0.01];

S = combinations(5, 0);
assert(isequal(cell2mat(S), zeros(1, 5)))
 
S = combinations(5, 1);
assert(isequal(cell2mat(S), eye(5)))

S = combinations(5, 5);
assert(isequal(cell2mat(S), ones(1, 5)))

S = combinations(3, 2);
S_test = ...
     [1     1     0;
      1     0     1;
      0     1     1];
assert(isequal(cell2mat(S), S_test))

S = combinations_lte(3, 3);
S_test = ...
     [0     0     0;
      1     0     0;
      0     1     0;
      0     0     1;
      1     1     0;
      1     0     1;
      0     1     1;
      1     1     1];
assert(isequal(cell2mat(S), S_test))

p = prob_Gammas(S, epsilon);
assert(size(p, 1) == size(S, 1))
assert(abs(p(1) - epsilon(1)^3) < 1e-10)
assert(all(abs(p(2:4) - epsilon(1)^2*epsilon(2)) < 1e-10))
assert(all(abs(p(5:7) - epsilon(1)*epsilon(2)^2) < 1e-10))
assert(abs(p(8) - epsilon(2)^3) < 1e-10)
assert(abs(sum(p) - 1) < 1e-15)

S = combinations_lte(3, 2);
assert(isequal(cell2mat(S), S_test(1:end-1,:)))

S = combinations_lte(5, 2);
p = prob_Gammas(S, epsilon);
assert(abs(1 - sum(p)) > 1e-10)

assert(isequal(cell2mat(combinations_lte(5, 0)), zeros(1, 5)))

S = combinations_gte_prob(3, 0, epsilon);
assert(isequal(cell2mat(S), zeros(1, 3)))

S = combinations_gte_prob(3, 0.97, epsilon);
assert(isequal(cell2mat(S), zeros(1, 3)))

S = combinations_gte_prob(3, 0.98, epsilon);
assert(isequal(cell2mat(S), [zeros(1, 3); eye(3)]))

S1 = combinations_gte_prob(5, 1, epsilon);
S2 = combinations_lte(5, 5);
assert(isequal(cell2mat(S1), cell2mat(S2)))

% Multi-variable sequences
S = combinations(6, 2);
S = cellfun(@(x) reshape(x, 2, []), S, 'UniformOutput', false);
epsilon = [0.99 0.95; 0.01 0.05];
p = prob_Gammas(S, epsilon);
assert(all(p(1) - 0.01*0.99*0.99*0.05*0.95*0.95 < 1e-10));

% Check for all S
for i = 1:numel(S)
    n_shocks = sum(S{i}, 2);
    n_no_shock = size(S{i},2) - n_shocks;
    prob_calc = prod(epsilon(2,:) .^ (n_shocks')) * ...
        prod(epsilon(1,:) .^ (n_no_shock'));
    assert(abs(prob_calc - p(i)) < 1e-12)
end

f = 3;
n_dist = 2;
m = 1;
d = 5;
S = combinations_lte(f*n_dist, m);
epsilon = [0.01; 0.01];
p = (ones(size(epsilon)) - (ones(size(epsilon)) - epsilon).^d)';
p = [ones(size(p))-p; p];
S = cellfun(@(x) reshape(x, n_dist, []), S, 'UniformOutput', false);
p_seq = prob_Gammas(S, p);
assert(abs(p_seq(1) - (p(1,1)^3*p(1,2)^3)) < 1e-10)
assert(all(abs(p_seq([2 4 6]) - (p(1,1)^2*p(2,1)*p(1,2)^3)) < 1e-10))
assert(all(abs(p_seq([3 5 7]) - (p(1,1)^3*p(1,2)^2*p(2,2))) < 1e-10))


%% Test most_probable_combinations.m

p = [0.95; 0.04; 0.01];
S = most_probable_sequences(p, 1, 0.99);
assert(isequal(S, [1; 2]))

p = [0.01 0.04 0.95];
S = most_probable_sequences(p, 1, 0.99);
assert(isequal(S, [3; 2]))

S = most_probable_sequences(p, 2, 0.99);
S_test = [   3     3;
             3     2;
             2     3;
             3     1];
for i=1:size(S,1)
    assert(isempty(setdiff(S(i,:), S_test(i,:))))
end

S = most_probable_sequences(p, 3, 0.95);
S_test = [   3     3     3;
             3     3     2;
             3     2     3];
for i=1:size(S,1)
    assert(isempty(setdiff(S(i,:), S_test(i,:))))
end
