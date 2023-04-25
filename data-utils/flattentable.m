function T = flattentable(T, label)
% T = flattentable(T)
% Create a table with one row containing all the elements of
% table t, by combining row and column names as columns of
% flattened table.
%
% Examples:
% >> names = ["Jill"; "Jack"]; age = [21; 32]; weight = [76; 81];
% >> T = table(age, weight, 'RowNames', names);
% >> flattentable(T)
% 
% ans =
% 
%   1x4 table
% 
%     T_Jill_age    T_Jill_weight    T_Jack_age    T_Jack_weight
%     __________    _____________    __________    _____________
% 
%         21             76              32             81  
% 
% >> flattentable(T, 'Data')
% 
% ans =
% 
%   1x4 table
% 
%     Data_Jill_age    Data_Jill_weight    Data_Jack_age    Data_Jack_weight
%     _____________    ________________    _____________    ________________
% 
%          21                 76                32                 81  
% 
%
    if nargin == 1
        % Use name of first argument if it has one
        label = inputname(1);
        if isempty(label)
            error("ValueError: Provide label for table.");
        end
    end

    % Extract data from table as array
    X = T.Variables;
    n = numel(X);
    [R, C] = ind2sub(size(X), 1:numel(X));
    row_names = T.Properties.RowNames';
    if isempty(row_names)
        row_names = string(1:size(X, 1));
    end
    row_names = row_names(R);
    var_names = T.Properties.VariableNames(C);
    col_names = join([repmat({char(label)}, 1, n); row_names; var_names], '_', 1);
    T = array2table(reshape(X,1,[]), 'VariableNames', col_names);

end