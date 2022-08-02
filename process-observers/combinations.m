function S = combinations(n, m)
% S = combinations(n, m) returns a cell array of
% integer vectors of width n representing all
% combinations of sequences with m values set to
% 1 and all other values 0.
%
% Example:
% >> S = combinations(3, 2);
% >> cell2mat(S)
% 
% ans =
% 
%   3×3 int16 matrix
% 
%    1   1   0
%    1   0   1
%    0   1   1
%      

    C = int16(nchoosek(1:n, m));
    S = repmat({int16(zeros(1, n))}, size(C, 1), 1);
    for i = 1:size(C, 1)
        S{i}(1, C(i,:)) = 1;
    end

end