function labels = matrix_element_labels(label,row_subs,col_subs,paren,lt)
% labels = matrix_element_labels(label,row_subs,col_subs,paren,lt)
% Makes a cell array matrix of labels based on the row and column
% labels provided.
%
% Example 1:
% >> matrix_element_labels('a', 1:2, 1:3, '')
% 
%   2×3 cell array
% 
%     {'a_1_1'}    {'a_1_2'}    {'a_1_3'}
%     {'a_2_1'}    {'a_2_2'}    {'a_2_3'}
% 
% Example 2:
% >> matrix_element_labels('x', {'a','b','c'}, {1,2}, '(t)')
% 
% ans =
% 
%   3×2 cell array
% 
%     {'x_a_1(t)'}    {'x_a_2(t)'}
%     {'x_b_1(t)'}    {'x_b_2(t)'}
%     {'x_c_1(t)'}    {'x_c_2(t)'}
% 
    if nargin == 4
        lt = false;
    end
    n_rows = numel(row_subs);
    n_cols = numel(col_subs);
    labels = cell(n_rows, n_cols);
    for r = 1:n_rows
        for c = 1:n_cols
            labels{r,c} = sprintf('%s_%s_%s%s', label, string(row_subs(r)), string(col_subs(c)), paren);
        end
    end
    if lt == true
        labels = string2latex(labels);
    end