function T = array2tablerow(X, label)
% T = array2tablerow(X, label)
% Create a table with one row containing all the elements of
% array X, using labels based on its name or the given label.
%
% Examples:
% >> X = 1:3;
% >> array2tablerow(X)
% 
% ans =
% 
%   1×3 table
% 
%     X_1    X_2    X_3
%     ___    ___    ___
% 
%      1      2      3 
% 
% >> array2tablerow(ones(3, 2), 'M')
% 
% ans =
% 
%   1×6 table
% 
%     M_1_1    M_2_1    M_3_1    M_1_2    M_2_2    M_3_2
%     _____    _____    _____    _____    _____    _____
% 
%       1        1        1        1        1        1  
% 

    if nargin == 1
        % Use name of first argument if it has one
        label = inputname(1);
        if isempty(label)
            error("ValueError: Provide label for array.");
        end
    end
    if isscalar(X)
        T = array2table(X, 'VariableNames', {label});
    else
        switch numel(size(X))

            case 2
                [R, C] = ind2sub(size(X), 1:numel(X));
                ind = [R' C'];
                n = size(ind, 1);
                col_names = cell(1, n);
                if ismember(1, size(X))
                    j = find(size(X) ~= 1);
                    for i = 1:n
                        col_names{i} = sprintf('%s_%d', label, ind(i, j));
                    end
                else
                    for i = 1:n
                        col_names{i} = sprintf('%s_%d_%d', label, ind(i,:));
                    end
                end

            case 3
                [R, C, D] = ind2sub(size(X), 1:numel(X));
                ind = [R' C' D'];
                n = size(ind, 1);
                col_names = cell(1, n);
                for i = 1:n
                    col_names{i} = sprintf('%s_%d_%d_%d', label, ind(i,:));
                end

        end
        T = array2table(reshape(X,1,[]), 'VariableNames', col_names);

    end
end