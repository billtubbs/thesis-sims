% Test downsample_evenly.m

%% Vector time series
X = (1:14)';
[Xmean, Xmin, Xmax] = downsample_evenly(X, 7);
assert(isequal(Xmean, [1.5 3.5 5.5 7.5 9.5 11.5 13.5]'))
assert(isequal(Xmin, [1 3 5 7 9 11 13]'))
assert(isequal(Xmax, [2 4 6 8 10 12 14]'))

rng(0)
X = ones(1000, 2);
[Xmean, Xmin, Xmax] = downsample_evenly(X(:, 1), 500);
assert(isequal(sum(Xmean), 500))
assert(isequal(Xmin, Xmean))
assert(isequal(Xmax, Xmean))


%% Multiple time series

A = reshape(1:24, 8, []);
[Xmean, Xmin, Xmax] = downsample_evenly(A, 4);
assert(isequal(Xmean, [...
    1.5000    9.5000   17.5000
    3.5000   11.5000   19.5000
    5.5000   13.5000   21.5000
    7.5000   15.5000   23.5000
]))
assert(isequal(Xmin, [...
     1     9    17
     3    11    19
     5    13    21
     7    15    23
]))
assert(isequal(Xmax, [...
     2    10    18
     4    12    20
     6    14    22
     8    16    24
]))

A = reshape(1:12, 4, []);

[Xmean, Xmin, Xmax] = downsample_evenly(A);
assert(isequal(Xmean, A))
assert(isequal(Xmin, A))
assert(isequal(Xmax, A))

[Xmean, Xmin, Xmax] = downsample_evenly(A, 3);
assert(isequal(Xmean, A))
assert(isequal(Xmin, A))
assert(isequal(Xmax, A))

[Xmean, Xmin, Xmax] = downsample_evenly(A, 2);
assert(isequal(Xmean, [...
    1.5000    5.5000    9.5000
    3.5000    7.5000   11.5000
]))
assert(isequal(Xmin, [...
     1     5     9
     3     7    11
]))
assert(isequal(Xmax, [...
     2     6    10
     4     8    12
]))

