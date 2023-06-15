% Test the following functions:
%
%  - prob_gamma.m
%  - prob_w_given_gamma.m
%  - prob_w.m
%  - prob_gamma_given_w.m
%  - n_filters.m
%  - seq_prob.m
%  - prob_gamma.m
%

clear all


%% Test probability funcs - one discrete binary variable

% Probability of discrete events
epsilon = [0.99; 0.01];

% Test prob_gamma.m
assert(prob_gamma(0, epsilon) == 0.99)
assert(prob_gamma(1, epsilon) == 0.01)
assert(isequal(prob_gamma([0 0 1], epsilon), [0.99 0.99 0.01]))

% Test prob_rk.m
assert(prob_rk(1, epsilon) == 0.99)
assert(prob_rk(2, epsilon) == 0.01)
assert(isequal(prob_rk([1 1 2], epsilon), [0.99 0.99 0.01]))

% Variance parameter
sigma1 = 0.01;
b = 100;
sigma_w = [sigma1; sigma1*b];

% Test prob_w_given_gamma.m
assert(round(prob_w_given_gamma(0, 0, sigma_w), 4) == 39.8942);
assert(round(prob_w_given_gamma(0, 1, sigma_w), 4) == 0.3989);
assert(round(prob_w_given_gamma(1, 0, sigma_w), 4) == 0.0);

% Test prob_w.m
assert(round(prob_w(0, epsilon, sigma_w), 4) == 39.4993)
assert(all(round(prob_w([0.1; 1; -10], epsilon, sigma_w), 4) == [0.0040; 0.0024; 0.0000]))
assert(abs(prob_w(0.01, epsilon, sigma_w) - ...
    (prob_w_given_gamma(0.01, 0, sigma_w) * epsilon(1) ...
    + prob_w_given_gamma(0.01, 1, sigma_w) * epsilon(2))) < 1e-12)

% Test prob_gamma_given_w.m
wk = 0; gamma_k = 0;
assert(round(prob_gamma_given_w(wk, gamma_k, epsilon, sigma_w), 4) == 0.9999)
wk = 1; gamma_k = 1;
assert(round(prob_gamma_given_w(wk, gamma_k, epsilon, sigma_w), 4) == 1)

% Tests with vector inputs
w = linspace(-3, 3, 101);
p1 = normpdf(w,0,sigma_w(1));
p2 = normpdf(w,0,sigma_w(2));
assert(max(abs(prob_w_given_gamma(w, 0, sigma_w) - p1)) < 1e-12);
assert(max(abs(prob_w_given_gamma(w, 1, sigma_w) - p2)) < 1e-12);
assert(max(abs(prob_w(w, epsilon, sigma_w) - (epsilon(2)*p2 + epsilon(1)*p1))) < 1e-12);

%% Test Markov model transition probabilities

% Binary variable
T = [0.95 0.05; 0.01 0.99];
gamma_km1 = [1 2 1 2]';
gamma_k =   [1 1 2 2]';
assert(isequal(prob_transitions(gamma_k, gamma_km1, T), ...
    [0.95 0.01 0.05 0.99]'))

% Discrete random variable with 3 modes
% m1 : Bull market
% m2 : Bear market
% m3 : Stagnant market
T = [0.975 0.02  0.005
     0.3   0.5   0.2;
     0.02  0.4   0.58];
r_km1 = [1 2 3 1 2 3 1 2 3]';
r_k   = [1 1 1 2 2 2 3 3 3]';
assert(isequal(prob_transitions(r_k, r_km1, T), ...
    [0.975 0.3 0.02 0.02 0.5 0.4 0.005 0.2 0.58]'))


%% Test probability funcs - two discrete binary variables

% Probability of discrete events
epsilon = [0.99 0.97; 0.01 0.03];

% Test prob_gamma.m
assert(isequal(prob_gamma([0; 0], epsilon), [0.99; 0.97]))
assert(isequal(prob_gamma([1; 1], epsilon), [0.01; 0.03]))
assert(isequal(prob_gamma([0 0 1; 1 0 0], epsilon), ...
    [0.99 0.99 0.01; 0.03 0.97 0.97]))

% Test prob_rk.m
assert(isequal(prob_rk([1; 1], epsilon), [0.99; 0.97]))
assert(isequal(prob_rk([2; 2], epsilon), [0.01; 0.03]))
assert(isequal(prob_rk([1 1 2; 2 1 1], epsilon), ...
    [0.99 0.99 0.01; 0.03 0.97 0.97]))

% Variance parameter
sigma1 = [0.01 0.005];
b = [100 100];
sigma_w = [sigma1; sigma1.*b];

% Test prob_w_given_gamma.m
assert(round(prob_w_given_gamma([0; 0], [0; 0], sigma_w), 1) == 3.1831e+03);
assert(round(prob_w_given_gamma([0; 0], [1; 1], sigma_w), 5) == 0.31831);
assert(round(prob_w_given_gamma([1; 1], [0; 0], sigma_w), 4) < 1e-15);


%% Test probability funcs - one discrete variable with 3 states

epsilon = [0.98 0.015 0.005]';

% Test prob_gamma.m
assert(prob_gamma(0, epsilon) == 0.98)
assert(prob_gamma(1, epsilon) == 0.015)
assert(prob_gamma(2, epsilon) == 0.005)
assert(isequal(prob_gamma([0 0 1 2], epsilon), [0.98 0.98 0.015 0.005]))

% Test prob_rk.m
assert(prob_rk(1, epsilon) == 0.98)
assert(prob_rk(2, epsilon) == 0.015)
assert(prob_rk(3, epsilon) == 0.005)
assert(isequal(prob_rk([1 1 2 3], epsilon), [0.98 0.98 0.015 0.005]))

% Variance parameter
sigma1 = 0.01;
b = [100 200];
sigma_w = [sigma1 sigma1*b]';

% Test prob_w_given_gamma.m
assert(round(prob_w_given_gamma(0, 0, sigma_w), 4) == 39.8942);
assert(round(prob_w_given_gamma(0, 1, sigma_w), 4) == 0.3989);
assert(round(prob_w_given_gamma(0, 2, sigma_w), 4) == 0.1995);
assert(round(prob_w_given_gamma(1, 0, sigma_w), 4) == 0.0);
assert(round(prob_w_given_gamma(1, 1, sigma_w), 4) == 0.2420);
assert(round(prob_w_given_gamma(1, 2, sigma_w), 4) == 0.1760);


%% Test probability funcs - two discrete variables with 3 states

epsilon = [0.98 0.015 0.005; 0.97 0.02 0.01]';

% Test prob_gamma.m
assert(isequal(prob_gamma([0; 0], epsilon), [0.98; 0.97]))
assert(isequal(prob_gamma([1; 1], epsilon), [0.015; 0.02]))
assert(isequal(prob_gamma([2; 2], epsilon), [0.005; 0.01]))
assert(isequal(prob_gamma([0 0 1 2; 2 1 0 0], epsilon), ...
    [0.98 0.98 0.015 0.005; 0.01 0.02 0.97 0.97]))

% Test prob_rk.m
assert(isequal(prob_rk([1; 1], epsilon), [0.98; 0.97]))
assert(isequal(prob_rk([2; 2], epsilon), [0.015; 0.02]))
assert(isequal(prob_rk([3; 3], epsilon), [0.005; 0.01]))
assert(isequal(prob_rk([1 1 2 3; 3 2 1 1], epsilon), ...
    [0.98 0.98 0.015 0.005; 0.01 0.02 0.97 0.97]))


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
for k = 0:m
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


% Example 3 - single disturbance used for report simulation
nw = 1;
epsilon = 0.01;
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
