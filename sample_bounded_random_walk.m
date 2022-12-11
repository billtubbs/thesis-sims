function p = sample_bounded_random_walk(sd_e, beta, alpha1, alpha2, N, ...
    tau, phi)
% p = sample_bounded_random_walk(sd_e, beta, alpha1, alpha2, N, ...
%     tau, phi)
% Simulate the Bounded Random Walk stochastic process proposed by
% J. Nicolau:
%
% Referemce:
%  - J. Nicolau, Stationary Processes That Look Like Random 
%    Walks - The Bounded Random Walk Process in Discrete and 
%    Continuous Time, Econometric Theory, 18, 2002, 99-118.
%

    % Regularization parameter (set to 0.5 by default)
    if nargin < 7
        phi = 0.5;
    end

    % Noise
    e = randn(N, 1);
    p = zeros(N, 1);

    % Set initial state
    xkm1 = tau;

    % Simulate
    for i = 1:N

        % Stochastic input
        alpha = sd_e * e(i);

        % Reversion bias
        bias = brw_reversion_bias(xkm1, alpha1, alpha2, beta, tau);

        % Regularization step (to avoid instability)
        if abs(bias) < 2 * abs(xkm1 - tau)
            x = xkm1 + bias + alpha;
        else
            x = tau + phi * (xkm1 - tau) + alpha;
        end

        % Bounded process output
        p(i) = x;
        xkm1 = x;

    end

end


