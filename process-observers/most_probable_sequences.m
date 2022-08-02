% function S = most_probable_sequences(p, n, p_cut)
%     assert(sum(p) == 1);
%     np = numel(p);
%     switch n
%         case 1
%             all_perms = (1:np)';
%         case 2
%             [X,Y] = ndgrid(1:np,1:np);
%             all_perms = [X(:) Y(:)];
%         case 3
%             [X,Y,Z] = ndgrid(1:np,1:np,1:np);
%             all_perms = [X(:) Y(:) Z(:)];
%         otherwise
%             error("n > 3 not implemented")
%     end
%     probs = prod(reshape(p(all_perms), [], n), 2);
%     [probs_sorted, order] = sort(probs, 'descend');
%     most_prob = cumsum(probs_sorted) <= p_cut;
%     S = all_perms(order(most_prob), :);
% end

% This answer from:
% https://www.mathworks.com/matlabcentral/answers/862870-efficient-algo
% rithm-to-compute-only-the-most-probable-sequences-of-a-random-variabl
% e-out-of-all-poss#answer_731520
%
% Requires uniqueperms from the file exchange:

function S = most_probable_sequences(p, n, p_cut)
    assert( abs(sum(p)-1) <= 1e-12);
    S =  recursor(p, n, p_cut);
    
end

function S = recursor(p, n, p_cut,~)
  if n == 1
    [p, S] = sort(p(:),'descend');   
    cut = cumsum(p) <= p_cut;
    S = S(cut);
    return
  end
  
   s = recursor(p, n-1, p_cut,[]);
   s = unique(sort(s, 2), 'rows');
   
   np = numel(p);
   ns = size(s,1);
   
   [I, J] = ndgrid(1:ns,1:np);
   
   s = unique(sort([s(I(:), :), J(:)], 2),'rows');
   
   s = num2cell(s,2);
   for i = 1:numel(s), s{i} = uniqueperms(s{i}); end
   s = cell2mat(s);
   
   sp = prod(p(s),2);
   [Probs, is] = sort(sp(:), 'descend');
  
   cut = cumsum(Probs) <= p_cut;
   
   S = s(is(cut),:);
   
end