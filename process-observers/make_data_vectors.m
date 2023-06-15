function vdata = make_data_vectors(vars, numtype)
% vdata = make_data_vectors(vars, numtype)
% Construct a set of vectors and associated meta data for
% the cell array of numeric variables passed to the function.
% This is useful for implementing S-functions in Simulink blocks 
% which require all working data to be stored as one or more 
% vectors.
%
% Example:
% >> vdata = make_data_vectors(1, [2 3; 4 5], {6, [7 8], 9})
% 
% vdata = 
% 
%   struct with fields:
% 
%      vecs: {[1]  [2 4 3 5]  [6 7 8 9]}
%     types: {'double'  'double'  {1×3 cell}}
%      dims: {[1 1]  [2 2]  {1×3 cell}}
% 
    if nargin < 2
        numtype = 'double';
    end
    n_vars = numel(vars);
    vecs = cell(1, n_vars);
    types = cell(1, n_vars);
    dims = cell(size(vars));
    n_els = cell(1, n_vars);
    for i = 1:n_vars
        arg = vars{i};
        switch class(arg)
            case numtype
                vecs{i} = reshape(arg, 1, []);
                types{i} = numtype;
                dims{i} = size(arg);
                n_els{i} = numel(arg);
            case 'cell'
                vd = make_data_vectors(arg, numtype);
                vecs{i} = cell2mat(vd.vecs);
                types{i} = reshape(vd.types, size(arg));
                dims{i} = reshape(vd.dims, size(arg));
                n_els{i} = numel(vecs{i});
            otherwise
                error("TypeError: invalid type.")
        end
    end
    vdata.vecs = vecs;
    vdata.types = types;
    vdata.dims = dims;
    vdata.n_els = n_els;
end