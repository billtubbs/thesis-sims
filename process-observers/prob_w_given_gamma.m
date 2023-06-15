function p = prob_w_given_gamma(wk, gamma_k, sigma_w)
    n = size(wk, 1);
    switch n
        case 1
            % Single variable
            p = normpdf(wk, 0, sigma_w(gamma_k+1));
        otherwise
            % Multi-variable
            s = sigma_w(sub2ind(size(sigma_w), gamma_k+1, (1:n)'));
            switch size(wk, 2)
                case 1
                    % Single set of values
                    p = mvnpdf(wk, [], diag(s.^2));
                otherwise
                    % Sequence of values
                    p = size(wk);
                    for i = 1:size(wk, 2)
                        p(:, i) = mvnpdf(wk(:, i), [], diag(s.^2));
                    end
            end
    end
return