function s2 = string2latex(s1)
% Converts a string or a cell array of strings to latex input
% format by adding the '$' character at the beginning and end
% of each string.
    switch class(s1)
        case {'char', 'string'}
            if ~(s1(1) == '$' && s1(end) == '$')
                s2 = strcat('$', s1, '$');
            else
                s2 = s1;
            end
        case 'cell'
            s2 = cellfun(@string2latex, s1, 'UniformOutput', false);
        otherwise
            error('object class not supported')
    end
end