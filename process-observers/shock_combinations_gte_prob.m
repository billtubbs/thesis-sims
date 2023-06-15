function S = shock_combinations_gte_prob(n, beta, epsilon)
% S = combinations_gte_prob(n, beta, epsilon) returns a cell 
% array of vectors of width n representing all combinations 
% of sequences which have a total probability of occurence
% greater than or equal to beta, given the probabilities in 
% epsilon.
%
% Example:
% >> S = shock_combinations_gte_prob(3, 0.99, [0.99; 0.01]);
% >> cell2mat(S)
% 
% ans =
% 
%   4×3 int16 matrix
% 
%    1   1   1
%    2   1   1
%    1   2   1
%    1   1   2
% 

    C = cell(n+1, 1);
    p_cum = 0;
    for i = 0:n
        C{i+1} = shock_combinations(n, i);
        p = prob_seq(C{i+1}, epsilon);
        p_cum = p_cum + sum(p);
        if p_cum >= beta
            break
        end
    end
    S = vertcat(C{:});
end