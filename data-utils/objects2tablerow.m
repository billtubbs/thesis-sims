function T = objects2tablerow(vars)
% T = objects2tablerow(vars) converts a containers.Map of 
% variables (vars) into a one-row table T using the keys
% of the Map as column headings. Handles objects of the 
% following classes:
%  - numeric (scalars and arrays)
%  - char
%  - string
%  - string array
%  - struct
%  - cell arrays
%
% Example:
% >> person.name = "Amy";
% >> person.age = int16(5);
% >> person.siblings = {'Peter', 'Fred'};
% >> data = [0.5377  1.8339 -2.2588];
% >> vars = containers.Map({'Person', 'Data'}, {person, data});
% >> T = objects2tablerow(vars)
% 
% T =
% 
%   1×7 table
% 
%     Data_1    Data_2    Data_3     Person_age    Person_name    Person_siblings_1    Person_siblings_2
%     ______    ______    _______    __________    ___________    _________________    _________________
% 
%     0.5377    1.8339    -2.2588        5            "Amy"           {'Peter'}            {'Fred'}     
% 

    var_names = vars.keys;
    n_vars = numel(var_names);
    sections = {};
    n = 0;
    for i = 1:n_vars

        if isnumeric(vars(var_names{i})) || islogical(vars(var_names{i}))
            n = n + 1;
            if numel(vars(var_names{i})) == 1
                sections{n} = table(vars(var_names{i}), 'VariableNames', var_names(i));
            else
                sections{n} = array2tablerow(vars(var_names{i}), var_names{i});
            end

        elseif ischar(vars(var_names{i}))
            n = n + 1;
            sections{n} = table({vars(var_names{i})}, 'VariableNames', var_names(i));

        elseif isstring(vars(var_names{i}))
            n = n + 1;
            if numel(vars(var_names{i})) == 1
                sections{n} = table(vars(var_names{i}), 'VariableNames', var_names(i));
            else
                sections{n} = array2tablerow(vars(var_names{i}), var_names{i});
            end

        elseif isstruct(vars(var_names{i}))
            n = n + 1;
            keys = fieldnames(vars(var_names{i}));
            n_items = numel(keys);
            for j = 1:n_items
                keys{j} = sprintf('%s_%s', var_names{i}, keys{j});
            end
            sections{n} = objects2tablerow(containers.Map(keys, struct2cell(vars(var_names{i}))));

        elseif iscell(vars(var_names{i}))
            n = n + 1;
            n_items = numel(vars(var_names{i}));
            keys = cell(1, n_items);
            for j = 1:n_items
                keys{j} = sprintf('%s_%d', var_names{i} , j);
            end
            sections{n} = objects2tablerow(containers.Map(keys, vars(var_names{i})));

        elseif istable(vars(var_names{i}))
            n = n + 1;
            n_items = numel(vars(var_names{i}));
            keys = cell(1, n_items);
            for j = 1:n_items
                keys{j} = sprintf('%s_%d', var_names{i} , j);
            end
            sections{n} = objects2tablerow(containers.Map(keys, vars(var_names{i})));

        else
            error('ValueError: Object class not supported')

        end
    end
    T = [sections{:}];
end