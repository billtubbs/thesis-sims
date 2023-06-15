function seq = indseq_from_drvs(Gamma, S)
% seq = indseq_from_drvs(Gamma)
% Converts a set of binary shock indicator sequences
% to a sequence of integer shock indicators where
% the integer values represent an index to the model
% corresponding to a unique combination of shocks.
% 
% Important notes:
%  - the values in seq are one-based (1, 2, ..., n) 
%    whereas those in Gamma are zero-based (0, 1).
%  - the convention used to represent combinations
%    is not what you might expect. See the following
%    scripts to understand how the integer values
%    relate to combinations of binary variables:
% 
% Related functions:
%  - shock_combinations.m
%  - shock_combinations_lte.m
%  - shock_combinations_gte_prob.m
%
% Arguments:
%  Gamma : (nT, nd) double
%      Set of nd binary shock indicator sequences 
%      arranged in columns. All values should be in
%      (0, 1).
%  S : (nj, 1) cell array (optional)
%      Cell array of possible combinations of 
%      system modes (1 = no shock, 2 = shock).
%
% Example 1:
% >> Gamma1 = [0 0 1 0 1]';
% >> Gamma2 = [0 1 0 0 1]';
% >> seq = indseq_from_drvs([Gamma1 Gamma2])
% 
% seq =
% 
%      1     3     2     1     4
%
% Example 2:
% >> Gamma1 = [0 0 1 0 1]';
% >> Gamma2 = [0 1 0 0 1]';
% >> Gamma3 = [0 0 0 1 0]';
% >> S = shock_combinations_lte(3, 2);
% >> seq = indseq_from_drvs([Gamma1 Gamma2 Gamma3], S)
% 
% seq =
% 
%      1     3     2     4     5
%

    nd = size(Gamma, 2);
    if nargin < 2
        % Assume all combinations are included
        S = shock_combinations_lte(nd, nd);
    end
    S = cell2mat(S);
    nT = size(Gamma, 1);
    seq = zeros(1, nT);
    for i = 1:nT
        seq(i) = find(ismember(S, Gamma(i, :) + 1, 'rows'));
    end

end