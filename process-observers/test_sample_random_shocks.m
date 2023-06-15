%% Unit tests for sample_random_shocks.m
clear all


%% single sequence - 1 distribution

p = 0.2;

% Get random values from rng
rng(0,'twister')
a = rand(6, 1) < p';
rn2 = randn(6, 1);
rn1 = randn(sum(a), 1);

rng(0,'twister')
x0 = sample_random_shocks(6, 1);
assert(all(size(x0) == [6 1]))
assert(sum(x0 ~= 0) == 6)

rng(0,'twister')
[x0, a0] = sample_random_shocks(6, 1);
assert(all(size(x0) == [6 1]))
assert(sum(a0 == 1) == 6)

rng(0,'twister')
[x1, a1] = sample_random_shocks(6, p);
assert(all(a1 == a))

rng(0,'twister')
sigma1 = 10;
[x2, a2] = sample_random_shocks(6, p, sigma1);
assert(isequal(x2, x1*10))
assert(all(a2 == a))


%% single sequence - 2 distributions

p = 0.2;
sigma1 = 10;
sigma2 = 2;

% Get random values from rng
rng(0,'twister')
a = rand(6, 1) < p';
rn2 = randn(6, 1);
rn1 = randn(sum(a), 1);

rng(0,'twister')
x1 = sample_random_shocks(6, p, sigma1, sigma2);
rng(0,'twister')
[x2, a2] = sample_random_shocks(6, p, sigma1, sigma2);
assert(isequal(x1, x2))
assert(all(size(x1) == [6 1]))
assert(all(a2 == a))
assert(all(x1(a) == sigma1*rn1(1:2)))
assert(all(x1(~a) == sigma2*rn2(~a)))


%% Multiple sequences - 1 distribution

rng(0,'twister')
[x0, a0] = sample_random_shocks(6, 1);

rng(0,'twister')
x = sample_random_shocks([3 2], 0.2);
x_test = [0 0; 0 0; x0(1) x0(2)];
assert(isequal(x, x_test))

rng(0,'twister')
x = sample_random_shocks([3 2], [0.2; 0.2]);
assert(isequal(x, x_test))

rng(0,'twister')
[x, a] = sample_random_shocks([3 2], 0.2, 10);
assert(isequal(x, x_test*10))
assert(all(size(a) == size(x)))
assert(isequal(x ~= 0, a ~= 0))

rng(0,'twister')
[x, a] = sample_random_shocks([3 2], [1; 1]);
assert(sum(sum(x ~= 0)) == 6)
assert(all(a == 1, [1 2]))

rng(0,'twister')
[x, a] = sample_random_shocks([3 2], [0; 1]);
assert(sum(sum(x ~= 0)) == 3)
assert(all(a(:,1) == 0))
assert(all(a(:,2) == 1))

rng(0,'twister')
s = [3 5];
[x, a] = sample_random_shocks(s, rand(5,1), rand(5,1));
assert(isequal(size(x), s))
assert(all(x(a == 0) == 0))

rng(0,'twister')
sigma1 = [10; 100];
[x, a] = sample_random_shocks([3 2], 0.8, sigma1);
x_test = [0 0; 0 x0(2)*100; x0(1)*10 x0(3)*100];
assert(isequal(x, x_test))
assert(all(x(a == 0) == 0))


%% Multiple sequences - 2 distributions

s = [10 2];
p = [0.6; 0.75];
sigma1 = [10; 100];
sigma2 = [2; 5];
[x, a] = sample_random_shocks(s, p, sigma1, sigma2);
assert(isequal(size(x), s))
assert(isequal(size(a), s))
assert(all(sum(x ~= 0), [1 2]))


%% Numerical statistical tests

s = [1e6 2];
p = [0.6; 0.75];
sigma1 = [10; 100];

rng(0,'twister')
[x, a] = sample_random_shocks(s, p, sigma1);
nz = x ~= 0;
assert(abs(mean(nz(:,1)) - p(1)) < 1e-3)
assert(abs(mean(nz(:,2)) - p(2)) < 1e-3)
assert((std(x(a(:, 1), 1)) / sigma1(1) - 1) < 1e-3)
assert((std(x(a(:, 2), 2)) / sigma1(2) - 1) < 1e-3)

sigma2 = [2; 5];
rng(0,'twister')
[x, a] = sample_random_shocks(s, p, sigma1, sigma2);
assert((std(x(a(:, 1), 1)) / sigma1(1) - 1) < 1e-3)
assert((std(x(a(:, 2), 2)) / sigma1(2) - 1) < 1e-3)
