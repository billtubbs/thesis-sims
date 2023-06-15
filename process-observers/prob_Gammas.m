function p = prob_Gammas(S, epsilon)
% p = prob_Gamma(S, epsilon) returns a vector of the 
% probabilities of each sequence in S, where S is a
% cell array containing sequences as row vectors, or
% matrices in the case of multi-variable sequences 
% arranged in rows.
%
% Example 1:
% >> S = combinations(3, 1);
% >> cell2mat(S)
% 
% ans =
% 
%      1     0     0
%      0     1     0
%      0     0     1
% 
% >> epsilon = [0.99; 0.01];
% >> p = prob_Gammas(S, epsilon)
% 
% p =
% 
%     0.0098
%     0.0098
%     0.0098
% 
% >> assert(all(p == 0.01*0.99*0.99))
% 
% Example 2 - multivariable:
% >> S = combinations(6, 2);
% >> S = cellfun(@(x) reshape(x, 2, []), S, 'UniformOutput', false);
% >> S{1}
% 
% ans =
% 
%      1     0     0
%      1     0     0
% 
% >> epsilon = [0.99 0.95; 0.01 0.05];
% >> p = prob_Gammas(S, epsilon)
% >> assert(p(1) == 0.01*0.99*0.99*0.05*0.95*0.95)

%TODO: Delete this with mkf_observer objects

    p = nan(numel(S), 1);
    for i = 1:numel(S)
        p(i) = prod(prob_gamma(S{i}, epsilon), [2 1]);
    end

end