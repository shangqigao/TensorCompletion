function z = ttt(x, y, dims_x, dims_y)
%TTT Tensor-times-tensor for htensor.
%
%  TTT(X,Y) computes the outer product for two htensor objects X and Y.
% 
%  TTT(X,Y,XDIMS,YDIMS) computes the contracted product of X and Y in the
%  modes specified by the row vectors XDIMS and YDIMS.
% 
%  TTT(X,Y,DIMS) is equivalent to calling TTT(X,Y,DIMS,DIMS).
%
%  Additionally, the dimension trees of X and Y must fulfill the following
%  requirements:
%
%  1. a) There is a node of X with a subtree comprising exactly the
%  dimensions in XDIMS, or in the complement of XDIMS. The same is
%  true for Y and YDIMS.
%
%  - OR -
%
%  1. b) All the dimensions of X are in XDIMS. There are two nodes in
%  Y, such that YDIMS is exactly the union of their subtrees or the
%  complements of their subtrees. The distance between these nodes is
%  exactly 3, i.e., there are two nodes between them. Or, vice versa,
%  the statement holds when swapping X, XDIMS with Y, YDIMS.
%
%  - AND -
%
%  2. The contracted parts of both dimension trees must have the same
%  structure.
%
%  More details on these requirements can be found in Section 3.5.2 of
%  [C. Tobler, PhD thesis, ETH Zurich, 2012].
%
%  In the case of complex tensors, the complex conjugate of X is used.
%
%  Examples
%  X = htenrandn([4 2 3]);
%  Y = htenrandn([3 4 2]);
%  Z = ttt(X,Y)     %<-- outer product of X and Y
%  Z = ttt(X,X,1:3) %<-- inner product of X with itself
%  Z = ttt(X,Y,[1 2 3],[2 3 1]) %<-- inner product of X & Y
%  Z = ttt(X,Y,[1 3],[2 1]) %<-- product of X & Y along specified dims
% 
%  See also TTT, HTENSOR/TTM.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt
%
% The structure of this function follows the function ttt from the Tensor
% Toolbox.


% Check number of arguments
if(nargin <= 1)
  error('Requires at least 2 arguments.');
elseif(nargin == 2)
  dims_x = [];
end

if(~exist('dims_y', 'var'))
  dims_y = dims_x;
end

% Check x, y
if(~isa(x, 'htensor') || ~isa(y, 'htensor'))
  error('X and Y must be of class htensor.');
end

% Check dims_x, dims_y
if( ~isindexvector(dims_x) || numel(unique(dims_x)) < numel(dims_x) || ...
    ~isindexvector(dims_y) || numel(unique(dims_y)) < numel(dims_y) )
  error(['DIMS_X and DIMS_Y must be vectors of positive integers,' ...
         ' without duplicate entries.']);
end

% Check tensor sizes
sz_x = size(x);
sz_y = size(y);
if(~isequal(sz_x(dims_x), sz_y(dims_y)))
  error('Tensor dimensions must agree.');
end

% Find out what case this is, and the subtrees in question:
compl_dims_x = setdiff(1:ndims(x), dims_x);
compl_dims_y = setdiff(1:ndims(y), dims_y);

roots_x = subtree_roots(x, dims_x);
compl_roots_x = subtree_roots(x, compl_dims_x);

roots_y = subtree_roots(y, dims_y);
compl_roots_y = subtree_roots(y, compl_dims_y);

if( any([numel(roots_x), numel(compl_roots_x)] == 1) && ...
    any([numel(roots_y), numel(compl_roots_y)] == 1) )
  
  % Case of one node each
  if(numel(roots_x) == 1)
    node_x = roots_x;
  else
    node_x = -compl_roots_x;
  end
  
  if(numel(roots_y) == 1)
    node_y = roots_y;
  else
    node_y = -compl_roots_y;
  end
  
  z = one_node(x, y, dims_x, dims_y, node_x, node_y);
  
elseif(numel(dims_x) == ndims(x) || numel(dims_y) == ndims(y))
  
  % Case of two nodes in one tensor is possible
  
  % For convenience, swap x and y if necessary
  if( numel(dims_y) ~= ndims(y) )
    tmp = x; x = y; y = tmp;
    tmp = dims_x; dims_x = dims_y; dims_y = tmp;
    tmp = roots_x; roots_x = roots_y; roots_y = tmp;
    tmp = compl_roots_x; compl_roots_x = compl_roots_y; compl_roots_y = tmp;
    x = conj(x);
    y = conj(y);
  end
  
  x_parent = x.parent;
  
  % Check whether there are two subtrees in x
  
  % Check for roots of dims_x:
  % The two nodes must be connected by one edge:
  % either one is the other's parent, or they are both children of
  % the root node.
  if(numel(roots_x) == 2)
    nodes = sort(x_parent(roots_x), 'ascend');
    
    if(x_parent(nodes(2)) == nodes(1) || ...
       x_parent(nodes(1)) == 1 && x_parent(nodes(2)) == 1)
      node_x_child = nodes(2);
    end
  end
  
  % Check for roots of complement of dims_x:
  % The two nodes must be connected by one edge:
  % either one is the other's parent, or they are both children of
  % the root node.
  if(~exist('node_x_child', 'var') && numel(compl_roots_x) == 2)
    nodes = sort(x_parent(compl_roots_x), 'ascend');
    
    if(x_parent(nodes(2)) == nodes(1) || ...
       x_parent(nodes(1)) == 1 && x_parent(nodes(2) == 1))
      node_x_child = nodes(2);
    end
  end
  
  if(~exist('node_x_child', 'var'))
    error(['Tensor-times-tensor cannot be calculated ' ...
	   'inside the hierarchical Tucker format, ' ...
	   'for this case.']);
  end
  
  % Case of two nodes in one tensor  
  z = two_nodes(x, y, dims_x, dims_y, node_x_child);
  
else
  error(['Tensor-times-tensor cannot be calculated ' ...
         'within the hierarchical Tucker format, ' ...
         'for the provided dimensions.']);
end


function z = one_node(x, y, dims_x, dims_y, node_x, node_y)
%
% Calculate ttt(x, y, dims_x, dims_y), given that
% 
%  node_x > 0: subtree of node_x contains all modes dims_x
%  node_x < 0: subtree of -node_x contains complement of dims_x
%
%  and that the same holds for node_y and dims_y.

% Change root such that right subtree contains all dims_x
if(node_x > 0)
  x = change_root(x, node_x, 'right');
else
  x = change_root(x, -node_x, 'left');
end

% Change root such that right subtree contains all dims_y
if(node_y > 0)
  y = change_root(y, node_y, 'right');
else
  y = change_root(y, -node_y, 'left');
end

% If all nodes or no nodes are chosen, node_x is +/-1 and change_root adds
% an additional node. All modes are shifted by one, and the new mode 1 has
% size 1. Adjust dims_x to this fact.
squeeze_left = false;
if(node_x == 1)
  dims_x = dims_x + 1;
  squeeze_left = true;
elseif(node_x == -1)
  dims_x = 1;
end

% If all nodes or no nodes are chosen, node_y is +/-1 and change_root adds
% an additional node. All modes are shifted by one, and the new mode 1 has
% size 1. Adjust dims_y to this fact.
squeeze_right = false;
if(node_y == 1)
  dims_y = dims_y + 1;
  squeeze_right = true;
elseif(node_y == -1)
  dims_y = 1;
end

M = elim_matrix(x, y, dims_x, dims_y);

% Combine x and y in one big tree, with some unused nodes

x.U = cellfun(@conj, x.U, 'UniformOutput', false);  
x.B = cellfun(@conj, x.B, 'UniformOutput', false);

offset_x = 1;
offset_y = x.nr_nodes + 1;

new_root_x = x.children(1, 1);
new_root_y = y.children(1, 1);

% Construct combined array children
x.children(x.children ~= 0) = x.children(x.children ~= 0) + offset_x;
y.children(y.children ~= 0) = y.children(y.children ~= 0) + offset_y;

children = [offset_x + new_root_x, offset_y + new_root_y; ...
            x.children; ...
            y.children];

% Construct combined cell arrays U and B
B{1} = M;
B(offset_x+(1:x.nr_nodes)) = x.B;
B(offset_y+(1:y.nr_nodes)) = y.B;

U(offset_x+(1:x.nr_nodes)) = x.U;
U(offset_y+(1:y.nr_nodes)) = y.U;

% Construct combined dims
dim2ind = [x.dim2ind + offset_x, y.dim2ind + offset_y];

% Eliminate unused nodes and adjust mode numbers
z = adjust_tree(children, dim2ind, U, B, x.is_orthog & y.is_orthog);

if(squeeze_left && squeeze_right)
  z = squeeze(z);
elseif(squeeze_right)
  z = squeeze(z, z.dims{z.children(1, 1)});
elseif(squeeze_left)
  z = squeeze(z, z.dims{z.children(1, 2)});
end

function z = two_nodes(x, y, dims_x, dims_y, node_x_child)
%
% Calculate ttt(x, y, dims_x, dims_y), given that dims_y = ndims(y) and
% dims_x belong to two subtrees, connected by one edge:
%
%   x
%    \
%    node_x_parent
%   /   \
%  x   node_x_child
%       /  \
%      x    x
%
% x: represent subtrees, two of these contain all dims_x

% Change root of x
x = change_root(x, node_x_child, 'right');

ii_left  = x.children(1, 1);
ii_right = x.children(1, 2);

% New x:
%
%         /  \
%      left  right
%      / \    / \
%     x   x  x   x
%

% Compute lr_left, lr_right such that:
% lr_left indicates which child of left contains some of dims_x,
% lr_right indicates which child of right contains some of dims_x.

% Go to bottom of left subtree of left
ii = ii_left;
while(~x.is_leaf(ii))
  ii = x.children(ii, 1);
end

if(any(dims_x == x.dims{ii}))
  lr_left = 1;
else
  lr_left = 2;
end

% Go to bottom of left subtree of right
ii = ii_right;
while(~x.is_leaf(ii))
  ii = x.children(ii, 1);
end

if(any(dims_x == x.dims{ii}))
  lr_right = 1;
else
  lr_right = 2;
end


% Split dims_x into the subtrees dims_x_left, dims_x_right and split dims_y
% accordingly
dims_x_left  = x.dims{x.children(ii_left , lr_left )};
dims_x_right = x.dims{x.children(ii_right, lr_right)};

dims_y_left = zeros(1, numel(dims_x_left));
for ii=1:numel(dims_x_left)
  dims_y_left(ii) = dims_y(dims_x == dims_x_left(ii));
end

dims_y_right = zeros(1, numel(dims_x_right));
for ii=1:numel(dims_x_right)
  dims_y_right(ii) = dims_y(dims_x == dims_x_right(ii));
end

% Check whether y can be split between dims_y_left and dims_y_right
roots_y_left  = subtree_roots(y, dims_y_left );
roots_y_right = subtree_roots(y, dims_y_right);

if(numel(roots_y_left) == 1)
  y = change_root(y, roots_y_left, 'left');
elseif(numel(roots_y_right) == 1)
  y = change_root(y, roots_y_right, 'right');
else
  error(['Tensor-times-tensor cannot be calculated ' ...
	 'inside the hierarchical Tucker format, ' ...
	 'for this case.']);
end

% Calculate elimination matrices:
%          
%         /    \
%    left        right       x
%   /   \       /     \
%     Mleft  Mright   
%         \   /              y
%
% (example for lr_left == 2, lr_right == 1)

x_left = change_root(x, ...
		     x.children(ii_left, lr_left), ...
		     'right');
M_left  = elim_matrix(x_left, y, dims_x_left, dims_y_left);

x_right = change_root(x, ...
		      x.children(ii_right, lr_right), ...
		      'right');
M_right = elim_matrix(x_right, y, dims_x_right, dims_y_right);

% Combine to form the overall elimination matrix

My = M_left * y.B{1} * M_right.';

B_ = ttm(x.B{ii_right}, {conj(x.B{1}), My}, [3, lr_right]);

M = ttt(x.B{ii_left}, B_, [3, lr_left], [3, lr_right]);

% Combine x and y in one big tree, with some unused nodes

x.U = cellfun(@conj, x.U, 'UniformOutput', false);  
x.B = cellfun(@conj, x.B, 'UniformOutput', false);

% Connect subtrees with root node, set root node matrix to M
x.B{1} = M;

x.children(1, 1) = x.children(ii_left , 3 - lr_left );
x.children(1, 2) = x.children(ii_right, 3 - lr_right);

% Eliminate unused nodes and adjust mode numbers
z = adjust_tree(x.children, x.dim2ind, x.U, x.B, x.is_orthog & y.is_orthog);


function roots = subtree_roots(x, dims)

v = false(x.nr_nodes, 1);

v(x.dim2ind(dims)) = true;

x_sibling = x.sibling;
x_parent = x.parent;

for ii=x.nr_nodes:-1:2
  
  if(v(ii) && v(x_sibling(ii)))
    v(x_parent(ii)) = true;
  end
    
end

roots = [];

for ii=x.nr_nodes:-1:2
  
  if(v(ii) && ~v(x_parent(ii)))
    roots(end+1) = ii;
  end
  
end

if(v(1))
  roots(end+1) = 1;
end


function M = elim_matrix(x, y, dims_x, dims_y)
%
% Calculate the elimination matrix from the right subtree of x, and
% whatever subtree of y it connects to.
%
% Note that dims_x, dims_y must each cover a subtree of x, y. This is not
% checked.

% Initialize mappings between nodes of x and y:
ix_leaves = x.dim2ind(dims_x);
iy_leaves = y.dim2ind(dims_y);

ix2iy = zeros(x.nr_nodes, 1);
iy2ix = zeros(y.nr_nodes, 1);

ix2iy(ix_leaves) = iy_leaves;
iy2ix(iy_leaves) = ix_leaves;

% Indices in right subtree of x
inds_x = htensor.subtree(x.children, x.children(1, 2));

M = cell(x.nr_nodes, 1);
x_is_leaf = x.is_leaf;
y_parent = y.parent;
y_is_left = y.is_left;

% Traverse the tree from the leaves upwards
for ix=inds_x(end:-1:1)
  
  if(x_is_leaf(ix))
    iy = ix2iy(ix);

    % M_t = U1_t' * U2_t
    M{ix} = full(x.U{ix}'*y.U{iy});
    
  else
    
    ix_left  = x.children(ix, 1);
    ix_right = x.children(ix, 2);
    
    iy_left  = ix2iy(ix_left);
    iy_right = ix2iy(ix_right);
    
    if(y_parent(iy_left) ~= y_parent(iy_right))
      error('Dimension trees of X and Y are incompatible.');      
    end
     
    iy = y_parent(iy_left);
    
    ix2iy(ix) = iy;
    iy2ix(iy) = ix;
    
    % M_t = B1_t' * (M_t1 kron M_t2) * B2_t
    B_ = ttm(x.B{ix}, { M{ix_left}, M{ix_right} }, [1 2], 't');
    
    szB_ = size(B_); szB_(end+1:3) = 1;
    szyB = size(y.B{iy}); szyB(end+1:3) = 1;
    sz_M = [szB_(3), szyB(3)];
    
    % Check whether ix_left is connected to iy_left or iy_right
    if(y_is_left(iy_left))
      M{ix} = ttt(B_, y.B{iy}, [1 2]);
    else
      M{ix} = ttt(B_, y.B{iy}, [1 2], [2 1]);
    end
    
    M{ix} = reshape(M{ix}, sz_M);
    
    % Save memory
    M{ix_left} = []; M{ix_right} = [];
    
  end      
end

M = M{x.children(1, 2)};


function x = adjust_tree(children, dim2ind, U, B, is_orthog)
%
% Eliminate all unused nodes, which are not in the subtree of node 1);
% rename modes to 1,...,d.
%

new2old = htensor.subtree(children, 1);

old2new = zeros(1, size(children, 1));
old2new(new2old) = 1:numel(new2old);

children = children(new2old, :);

is_leaf = (children(:, 1) == 0);

children(~is_leaf, :) = old2new(children(~is_leaf, :));

U = U(new2old);
B = B(new2old);

dim2ind = old2new(dim2ind);
dim2ind = dim2ind(dim2ind ~= 0);

x = htensor(children, dim2ind, U, B, is_orthog);
