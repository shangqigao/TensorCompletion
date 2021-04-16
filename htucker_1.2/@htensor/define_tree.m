function [children, dim2ind] = define_tree(dims, tree_type)
%DEFINE_TREE Define dimension tree.
%
%   [CHILDREN, DIM2IND] = DEFINE_TREE(DIMS) returns two arrays defining a
%   balanced dimension tree for an HTD.
%   DIMS specifies the order of the dimensions used when constructing the
%   tree. For an order-D tensor, DIMS = 1:D or any permutation thereof.
%   CHILDREN is an array of size nr_nodes x 2, containing in each row
%   the indices of the child nodes of the corresponding tree node. For leaf
%   nodes, both indices are zero. The index of a node is always
%   larger than the index of its parent node.
%   DIM2IND is a 1 x D vector where DIM2IND(i) is the index of the (leaf)
%   node that represents the i-th dimension.
%
%   DEFINE_TREE(DIMS, TREE_TYPE) allows to choose other, unbalanced trees.
%   Currently, the following options for TREE_TYPE are supported:
%
%   TREE_TYPE = 'first_separate': The left child of the root node is a leaf
%               (corresponding to the leading dimension in DIMS) and the
%               right subtree is balanced.
%   TREE_TYPE = 'first_pair_separate': The left subtree of the root node 
%               has two leaf nodes (corresponding to the leading two
%               dimensions in DIMS) and the right subtree is balanced.
%   TREE_TYPE = 'TT': Degenerate tree, corresponding to the Tensor Train
%               format.
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin < 1)
  error('Requires at least 1 argument.')
elseif( ~isindexvector(dims) || numel(unique(dims)) ~= numel(dims) )
  error('First argument must be an index vector without duplicate entries.')
end

if(nargin == 1)
  tree_type = '';
elseif(~ischar(tree_type))
  error('TREE_TYPE must be a character array.')
end

d = length(dims);
if (d == 0),
  error('Zero-dimensional tensors are not supported.');
end

% Check input dimensions
if( any(sort(dims) ~= 1:d) )
  error('All dimensions from 1 to d must be used in DIMS.')
end

children = zeros(2*d-1, 2);
dims = {dims};

nr_nodes = 1;
ii = 1;

while(ii <= nr_nodes)

  if(length(dims{ii}) == 1)
    children(ii, :) = [0, 0];
  else
    ii_left  = nr_nodes + 1;
    ii_right = nr_nodes + 2;
    nr_nodes = nr_nodes + 2;
    
    children(ii, :) = [ii_left, ii_right];
    
    if(nargin == 2 && strcmp(tree_type, 'first_separate') && ii==1)
      dims{ii_left } = dims{ii}(1);
      dims{ii_right} = dims{ii}(2:end);      
    elseif(nargin == 2 && strcmp(tree_type, 'first_pair_separate') && ii==1 && d>2)
      dims{ii_left } = dims{ii}(1:2);
      dims{ii_right} = dims{ii}(3:end);      
    elseif(nargin == 2 && strcmp(tree_type, 'TT'))
      dims{ii_left } = dims{ii}(1);
      dims{ii_right} = dims{ii}(2:end);      
    else
      dims{ii_left } = dims{ii}(1:floor(end/2));
      dims{ii_right} = dims{ii}(floor(end/2)+1:end);
    end
  end
  
  ii = ii + 1;
end

ind_leaves = find(children(:, 1) == 0);

dim2ind(cell2mat(dims(ind_leaves))) = ind_leaves;
