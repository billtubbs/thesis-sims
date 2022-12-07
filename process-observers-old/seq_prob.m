function p_seq = seq_prob(m, n, p)
% pr = seq_prob(m, n, p) returns the probability
% of no more than m shocks in a sequence of n 
% periods where the probability of a shock per
% period is p.

    % Cumulative binomial probability density
    p_seq = binocdf(m, n, 1-p);

end
