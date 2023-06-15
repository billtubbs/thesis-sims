% Test calc
% Method of calculating multi-model observer state estimates and
% error covariances matrices using 3-d arrays.  The code in
% MKFObserver.update() method is based on this
%

clear all

% Input test data

% State estimates of each filter
Xkf_est = cat(3, [1; 2], [3; 4], [5; 6]);

% Error covariance matrices of each filter estimate
Pkf_est = reshape(1:12, 2, 2, 3);

% Conditional probabilities of filter hypotheses
p_seq_g_Yk = [0.2; 0.5; 0.3];

% Convert probabilities into weights in 3rd dimension
weights = reshape(p_seq_g_Yk, 1, 1, []);

% Merged state estimate
xkp1_est = sum(Xkf_est .* weights, 3);
assert(isequal(xkp1_est, ([1 2] .* 0.2 + [3 4] .* 0.5 + [5 6] .* 0.3)'))
assert(isequal(round(xkp1_est, 4), [3.2 4.2]'))

% Deviations from merged estimates
Xkf_devs = xkp1_est - Xkf_est;
assert(isequal(round(Xkf_devs, 4), cat(3, [2.2; 2.2], [0.2; 0.2], [-1.8; -1.8])))

% Calculation of merged error covariance matrix
P = (Pkf_est + pagemtimes(Xkf_devs, pagetranspose(Xkf_devs))) .* weights;

test_P = cat(3, ...
    [1.1680    1.5680;
     1.3680    1.7680], ...
    [2.5200    3.5200;
     3.0200    4.0200], ...
    [3.6720    4.2720;
     3.9720    4.5720]);

assert(isequal(round(P, 4), test_P))
