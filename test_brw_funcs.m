% Test sample_bounded_random_walk.m and brw_reversion_bias.m

clear all
rng(0)

% Example from paper
beta = -15;  % notation change beta = k from Nicolau paper
alpha1 = 3;
alpha2 = 3;
tau = 100;

x = linspace(94, 106, 11);
% Values from Python code Bounded-random-walk-demo.ipynb
a_test = [
      2.00855369e+01  5.48811636e-01  1.49955768e-02  4.09734751e-04 ...
      1.11871265e-05  0.00000000e+00 -1.11871265e-05 -4.09734751e-04 ...
      -1.49955768e-02 -5.48811636e-01 -2.00855369e+01
];
a = brw_reversion_bias(x, alpha1, alpha2, beta, tau);
assert(max(abs(a - a_test)) < 1e-6)


% Example from paper
beta = -15;
alpha1 = 3;
alpha2 = 3;
tau = 100;
sd_e = 0.4;
N = 2000;

p = sample_bounded_random_walk(sd_e, beta, alpha1, alpha2, N, tau);
assert(all(p < 106) & all(p > 94))

assert(isequal(round(p(1:10), 4), [ ...
    100.2151  100.9486  100.0451  100.3899  100.5175 ...
    99.9944   99.8209   99.9580  101.3893  102.4971
]'))
