function vars = unpack_data_vectors(vdata)
% vars = unpack_data_vectors(vdata)
% Unpack a set of variables using the data in vdata.
% vdata is a struct created by the function 
% make_data_vectors for storing numerical data as
% a set of vectors, which is useful for implementing
% S-functions in Simulink blocks.
%
% Example:
% >> vdata = make_data_vectors(1, [2 3; 4 5], {6, [7 8], 9});
% >> [a, b, c] = unpack_data_vectors(vdata)
% 
% a =
% 
%      1
% 
% 
% b =
% 
%      2     3
%      4     5
% 
% 
% c =
% 
%   1×3 cell array
% 
%     {[6]}    {1×2 double}    {[9]}
% 
    n_vars = numel(vdata.vecs);
    vars = cell(size(vdata.dims));
    for i = 1:n_vars
        vec = vdata.vecs{i};
        type = vdata.types{i};
        dim = vdata.dims{i};
        if iscell(type)
            vd.vecs = mat2cell(vec, 1, cellfun(@numel_recursive, dim));
            vd.types = type;
            vd.dims = dim;
            vars{i} = unpack_data_vectors(vd);
        else
            vars{i} = reshape(vec, dim);
        end
    end