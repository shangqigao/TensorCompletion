function check_htensor(x)
%CHECK_HTENSOR Check internal consistency of htensor.
%
%   CHECK_HTENSOR(X) returns an error message if there is an inconsistency
%   in the fields of the htensor X.
%
%   This function is used by the constructor HTENSOR and can be used for
%   debugging purposes.
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Check consistency of children array.
if(~isnumeric(x.children) || size(x.children, 2) ~= 2 || ...
    any(floor(x.children(:)) ~= ceil(x.children(:))) || ...
    any(x.children(:) > size(x.children, 1)) || ...
    any(x.children(:) < 0) )
  error(['CHILDREN must be an N x 2 array of integers between' ...
	 ' 0 and N, defining a binary tree.']);
end

if(~isequal(sort(htensor.subtree(x.children, 1)), 1:x.nr_nodes))
  error(['Some indices of CHILDREN are not part of the tree, or' ...
	 ' the root node does not have index 1.']);
end

% Check that children indices are higher than parent indices and that all
% nodes have either zero or two children.
ind1 = find(x.children(:, 1) ~= 0);
ind2 = find(x.children(:, 2) ~= 0);
if(~isequal(ind1, ind2))
  error('Each node should have either zero or two children.');
end
tmp = x.children - (1:x.nr_nodes)'*ones(1, 2);
if(tmp(ind1, :) <= 0)
  error('Child nodes must have higher indices than their parent node.');
end

% Check consistency of dim2ind
if( ~isindexvector(x.dim2ind) || size(x.dim2ind, 1) ~= 1 )
  error('DIM2IND must be an index row vector.');
end

d = numel(x.dim2ind);

if( size(x.children, 1) ~= 2*d - 1 )
  error('Inconsistent number of nodes for CHILDREN and DIM2IND.');
end

if( any(x.children(x.dim2ind, 1) ~= 0) )
  error('Indices in DIM2IND associated with non-leaf nodes');
end

% Check data types of U and B
if(~isa(x.U, 'cell') || ...
   ~isa(x.B, 'cell') || ...
   ~all(cellfun(@(x)(isa(x, 'numeric')), x.U)) || ...
   ~all(cellfun(@(x)(isa(x, 'numeric')), x.B)) )
  error('U and B must be cell arrays of multidimensional arrays.');
end

if( size(x.U, 1) ~= 1 || size(x.B, 1) ~= 1 )
  error('U and B must be cell arrays of size 1 x nr_nodes.');
end

% Check ndims of U and B
if( ~all(cellfun('ndims', x.U) == 2) )
  error('Each U{ii} must be a matrix.');
end

if( ~all(ismember(cellfun('ndims', x.B), [2 3])) )
  error('Each B{ii} must be 2- or 3-dimensional array.');
end

% Check sizes of U and B
x_parent = x.parent;
x_is_left = x.is_left;
x_is_leaf = x.is_leaf;
for ii=2:size(x.children, 1)
  ii_par = x_parent(ii);
  if(x_is_left(ii))
    left_right = 1;
  else
    left_right = 2;
  end
  
  if(x_is_leaf(ii))
    if(size(x.U{ii}, 2) ~= size(x.B{ii_par}, left_right))
      error('Inconsistent sizes of U and B.');
    end
    if(size(x.U{ii}, 2) == 0)
      error('All hierarchical ranks must be at least 1.');
    end
  else
    if(size(x.B{ii}, 3) ~= size(x.B{ii_par}, left_right))
      error('Inconsistent sizes of U and B.');
    end      
    % Check that no rank is 0.
    if(any(size(x.B{ii}) == 0))
      error('All hierarchical ranks must be at least 1.');
    end
  end
end

% Check that root node has rank 1.
if(size(x.B{1}, 3) ~= 1)
  error('Root node must have hierarchical rank 1.');
end

% Check is_orthog
if( ~isa(x.is_orthog, 'logical') || ~isscalar(x.is_orthog) )
  error('IS_ORTHOG must be a logical scalar.');
end
