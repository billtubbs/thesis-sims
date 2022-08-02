function t = drop_empty_cols(t)
% t = drop_empty_cols(t)
% Use this to convert a table before writing it to a csv
% file using MATLAB's writetable function if the table may 
% contain empty arrays in one or more columns.
% This overcomes a bug in writetable where columns with empty
% arrays are missed from the CSV file even though the column
% names remain, resulting in mis-aligned data.
%
% Example
% >> U = randn(5, 0);  % empty array
% >> P = randn(5, 1);
% >> data = table(U, P);
% >> writetable(drop_empty_cols(data), 'mydata.csv')
% >> % This is equivalent to:
% >> writetable(data(:, 'P'), 'mydata.csv')
%
    varnames = t.Properties.VariableNames;
    n_cols = size(t, 2);
    cols_to_keep = boolean(ones(1, n_cols));
    for i = 1:n_cols
        name = varnames{i};
        if size(t{:, name}, 2) == 0
            cols_to_keep(i) = false;  % drop column
        end
    end
    t = t(:, varnames(cols_to_keep));
end