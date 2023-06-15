function p = prob_seq(S, epsilon)
% p = prob_seq(S, epsilon) returns a vector of the 
% probabilities of each sequence in S, where S is a
% cell array containing sequences as row vectors, or
% matrices in the case of multi-variable sequences 
% arranged in rows.
%
% Example 1:
% >> S = shock_combinations(3, 1);
% >> cell2mat(S)
% 
% ans =
% 
%   3×3 int16 matrix
% 
%    2   1   1
%    1   2   1
%    1   1   2
% 
% >> epsilon = [0.99; 0.01];
% >> p = prob_seq(S, epsilon)
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
% >> S = shock_combinations(6, 2);
% >> S = cellfun(@(x) reshape(x, 2, []), S, 'UniformOutput', false);
% >> S{1}
% 
% ans =
% 
%   2×3 int16 matrix
% 
%    2   1   1
%    2   1   1
% 
% >> epsilon = [0.99 0.95; 0.01 0.05];
% >> p = prob_seq(S, epsilon);
% >> assert(p(1) == 0.01*0.99*0.99*0.05*0.95*0.95)
%

    p = nan(numel(S), 1);
    for i = 1:numel(S)
        p(i) = prod(prob_rk(S{i}, epsilon), [2 1]);
    end

end