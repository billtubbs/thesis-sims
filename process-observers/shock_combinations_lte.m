function S = shock_combinations_lte(n, m)
% S = shock_combinations_lte(n, m) returns a cell array of
% integer vectors of width n representing all
% combinations of sequences with no more than m 
% values set to 2 and all other values 1.
%
% Example:
% S = shock_combinations_lte(3, 2);
% cell2mat(S)
% 
% ans =
% 
%   7×3 int16 matrix
% 
%    1   1   1
%    2   1   1
%    1   2   1
%    1   1   2
%    2   2   1
%    2   1   2
%    1   2   2
% 
    arguments
        n {mustBeInteger, mustBeGreaterThan(n, 0)}
        m {mustBeInteger}
    end

    C = cell(m+1, 1);
    for i = 0:m
        C{i+1} = shock_combinations(n, i);
    end
    S = vertcat(C{:});
end