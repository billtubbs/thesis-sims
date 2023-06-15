function S = shock_combinations(n, m)
% S = shock_combinations(n, m) returns a cell array of
% integer vectors of width n representing all
% combinations of sequences with m values set to
% 2 and all other values 1.
%
% Example:
% >> S = shock_combinations(3, 2);
% >> cell2mat(S)
% 
% ans =
% 
%   3×3 int16 matrix
% 
%    2   2   1
%    2   1   2
%    1   2   2
%    
    arguments
        n {mustBeInteger, mustBeGreaterThan(n, 0)}
        m {mustBeInteger, mustBeGreaterThanOrEqual(m, 0), ...
            mustBeLessThanOrEqual(m, n)}
    end

    if n > 1
        C = int16(nchoosek(1:n, m));
    else
        C = ones(1, m);
    end
    S = repmat({int16(ones(1, n))}, size(C, 1), 1);
    if ~isempty(C)
        for i = 1:size(C, 1)
            S{i}(1, C(i,:)) = 2;
        end
    end

end