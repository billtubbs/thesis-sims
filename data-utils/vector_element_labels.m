
function labels = vector_element_labels(label,paren,n,lt)
% labels = vector_element_labels(label,paren,n,lt)
% Creates a cell array of labels of length n suitable
% for labelling the columns of a matrix.
%
% Example
% >> vector_element_labels('x', '', 3)
% 
% ans =
% 
%   1×3 cell array
% 
%     {'x_1'}    {'x_2'}    {'x_3'}
% 
% >> vector_element_labels('\hat{x}', '(t)', 3, true)
% 
% ans =
% 
%   1×3 cell array
% 
%     {'$\hat{x}_1(t)$'}    {'$\hat{x}_2(t)$'}    {'$\hat{x}_3(t)$'}
% 
    if nargin == 3
        lt = false;
    end
    if n == 1
        labels = {strcat(label, paren)};
    else
        labels = cell(1, n);
        for i = 1:n
            labels{i} = sprintf('%s_%d%s', label, i, paren);
        end
    end
    if lt == true
        labels = string2latex(labels);
    end