function s2 = escape_latex_chars(s1)
% Converts a plain text string so that it can be used
% with the Latex interpreter in plot labels, legends, etc.
% by hiding superscript, subscript, backslash characters
% from the interpreter so they are not interpreted in the
% usual way.
%
% Example:
% >> escape_latex_chars("max_speed")
% 
% ans = 
% 
%     "max\_speed"
%
% Based on code from this answer:
% - https://www.mathworks.com/matlabcentral/answers/94238-how-can-i-place-the-_-or-characters-in-a-text-command#answer_737913
%
    switch class(s1)
        case {'char', 'string'}
            s2 = regexprep(s1, '[\\\^\_%]', '\\$0');
        case 'cell'
            s2 = cellfun(@escape_latex_chars, s1, 'UniformOutput', false);
        otherwise
            error('object class not supported')
    end
end