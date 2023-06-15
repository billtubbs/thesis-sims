function n = numel_recursive(dims)
% n = numel_recursive(dims)
% Recursively calculates the number of elements of a cell
% array described by the dimensions in dims. Dims may be
% a double describing the dimensions of one array or a cell
% array of dimensions, or a nested cell array of dimensions.
%
% This function is used by the function unpack_data_vectors
% to determine the correct size of the arrays to be unpacked.
%
% Example:
% >> numel_recursive({[2 1], [3 1], {[1 2], [1 3]}})
% 
% ans =
% 
%     10
% 
    if iscell(dims)
        n = sum(cellfun(@numel_recursive, dims));
    elseif isnumeric(dims)
        n = prod(dims);
    end
end