function x = change_root(x, ind, lr_subtree)
%CHANGE_ROOT Change root of the dimension tree.
%
%   Y = CHANGE_ROOT(X, IND) for an integer IND between 1 and 2*ndims(X)-1
%   changes the root of the dimension tree of the htucker X as follows: The
%   node IND becomes the right child of the root in Y (hence, its index
%   changes from IND to 3) and the rest of the tree is adjusted
%   accordingly.
%
%   Y = CHANGE_ROOT(X, IND, LR_SUBTREE) additionally allows to choose
%   whether IND becomes the right child (LR_SUBTREE = 'right', default) or
%   the left child (LR_SUBTREE = 'left') of the dimension tree.
%
%   Note that IND = 1 (root node) results in an additional level (and in a
%   singleton dimension):
%
%               LR_SUBTREE = 'right'   LR_SUBTREE = 'right'
%                        1                      1
%                      /   \                  /   \
%                    1    original       original   1
%                         dim tree       dim tree
%
%   Example:
%   x = htensor([4 6 2 3]);
%   y = change_root(x, 6);
%   size(y) %<-- returns [4 6 2 3]
%   z = change_root(x, 1);
%   size(z) %<-- returns [1 4 6 2 3]
%
%   Resulting dimension trees:
%   x:           y:   / \     z:   /  \
%       / \         / \  3        1   / \
%     /\   /\      /\  4            /\   /\ 
%    1  2 3  4    1  2             2  3 4  5
%
%   See also HTENSOR, CHANGE_DIMTREE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

if(nargin == 1)
  error('Requires at least 2 arguments.')
end

if(~isindexvector(ind) || ~isscalar(ind) || (ind > 2*ndims(x)-1))
  error('Second argument IND must be an integer between 1 and 2*ndims(x)-1.')
end

if(nargin == 2)
  lr_subtree = 'right';
elseif(~ischar(lr_subtree))
  error('LR_SUBTREE must be a string.')
end


if(ind ~= 1)
  
  children = x.children;
  p = x.parent;
  nr_nodes = x.nr_nodes;
  ind_par = p(ind);
  
  % Easy case: only need to switch children of the root node, if necessary
  if(ind_par == 1)
    if(strcmp(lr_subtree, 'left'))
      if(children(1, 2) == ind)
        children(1, [1 2]) = children(1, [2 1]);
      end
    else
      if(children(1, 1) == ind)
        children(1, [1 2]) = children(1, [2 1]);
      end
    end
  else
    
    % Change root node to connect ind and ind_par
    root_children = children(1, :);
    if(strcmp(lr_subtree, 'left'))
      children(1, 1) = ind;
      children(1, 2) = ind_par;
    else
      children(1, 2) = ind;
      children(1, 1) = ind_par;
    end
    
    % Change parent-child direction between ancestors of ind
    ii = ind;
    ii_par = ind_par;
    
    while(ii_par ~= 1)
      
      left_right = find(children(ii_par, :) == ii);
      ii = ii_par;
      ii_par = p(ii);
      
      children(ii, left_right) = ii_par;
    
    end
    
    % Reconnect children of the original root node
    children(ii, left_right) = root_children(root_children ~= ii);
    
  end
  
  % Rearrange indices, such that 1:nr_nodes descends through the tree
  new2old = ones(1, nr_nodes);
  ii_read = 1;
  ii_write = 2;
  
  while(ii_write < nr_nodes)
    
    if( all(children(new2old(ii_read), :) ~= [0 0]) )
      new2old(ii_write:ii_write+1) = children(new2old(ii_read), :);
      ii_write = ii_write + 2;
    end
    
    ii_read = ii_read+1;
  end
  
  old2new(new2old) = 1:nr_nodes;
  
  % Construct children and dim2ind of new htensor
  children = children(new2old, :);
  for ii=1:size(children, 1)
    if( any(children(ii, :) ~= 0) )
      children(ii, :) = old2new(children(ii, :));
    end
  end
      
  dim2ind = old2new(x.dim2ind);
  
  x = change_dimtree(x, children, dim2ind);
    
else % if ind == 1 (root node)                1
     % create additional level above root:   /  \
     %                                      1  orig_tree
     
  children = x.children;
  children(children ~= 0) = children(children ~= 0) + 2;
  
  if(strcmp(lr_subtree, 'left'))
    children = [3 2; 0 0; children];
  else  
    children = [2 3; 0 0; children];
  end
  
  dim2ind = [2, x.dim2ind + 2];
  
  U((1:x.nr_nodes)+2) = x.U;
  U{2} = 1;
  
  B((1:x.nr_nodes)+2) = x.B;
  B{1} = 1;
  
  x = htensor(children, dim2ind, U, B, false);
  
end
