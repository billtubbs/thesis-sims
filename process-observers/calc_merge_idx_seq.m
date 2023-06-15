function merge_idx_seq = calc_merge_idx_seq(seq)
% merge_idx_seq = calc_merge_idx_seq(seq)
% Calculates a sequence of column vectors which indicate
% which sequences should be merged at the end of each detetion
% interval in the sequence fusion algorithm described by
% Robertson et al. (1998), see p. 265:
%  "those sequences (δ1, δ2, ..., δk} whose last h, say, terms 
%   are the same can be merged into a single sequence (i.e. we
%   do not discriminate among those sequences which have the
%   last h terms in common)."
%
% Returns a matrix of int16 index values. Each column i shows
% the sequences which are unique over all detection intervals
% except i.  See example below. When i = 1, sequences 1 and 2
% are the same at intervals 2 and 3 (both contain [0 0]). 
% Sequences 3 and 4 are unique (they contain [1 0] and [0 1]
% at intervals 2 and 3).
%
% Example:
% >> seq = {[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]};
% >> merge_idx_seq = calc_merge_idx_seq(seq)
% 
% merge_idx_seq =
% 
%   4×3 int16 matrix
% 
%    1   1   1
%    1   2   2
%    2   1   3
%    3   3   1
%
    sz = size(seq{1});  % assume all the same size
    assert(sz(1) == 1);
    seq = cell2mat(seq);
    merge_idx_seq = int16(zeros(size(seq)));
    for i = 1:sz(2)
        seq_fuse = [seq(:,i+1:end) seq(:,1:i-1)];
        [~,~,merge_idx] = unique(seq_fuse,'rows','stable');
        merge_idx_seq(:, i) = merge_idx;
    end
end