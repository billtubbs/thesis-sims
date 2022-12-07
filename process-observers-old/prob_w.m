function p = prob_w(wk, epsilon, sigma_w)
    assert(all((epsilon >= 0) & (epsilon <= 1)))
    switch numel(wk)
        case 1
            p = sum(epsilon .* normpdf(wk, 0, sigma_w), 1);
        otherwise
            p = zeros(size(wk));
            for i = 1:numel(wk)
                p(i) = sum(epsilon .* normpdf(wk(i), 0, sigma_w), 1);
            end
    end
return