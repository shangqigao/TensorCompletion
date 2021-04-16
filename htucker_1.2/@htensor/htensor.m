classdef htensor
%HTENSOR - Hierarchical Tucker tensor.
%
%   A MATLAB class for representing a tensor in Hierarchical Tucker
%   Decomposition (HTD), allowing to construct, manipulate and operate with
%   such tensors.
%
%   See also HTENSOR/HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

properties( SetAccess = private, GetAccess = public )

children   % Structure of the dimension tree.
dim2ind    % Index of leaf node for each dimension.
U          % Matrices in the leaf nodes.
B          % Transfer tensors in non-leaf nodes.
is_orthog  % Indicates whether the HTD is orthogonalized.

end

% Dependent properties
properties( Dependent = true, SetAccess = private, GetAccess = public )

parent     % Parent of each node.
nr_nodes   % Number of nodes in tree.
is_leaf    % Indicates for each node whether it is a leaf.
lvl        % Level (distance from root) of each node.
dims       % Dimensions represented by each node.
bin_dims   % Dimensions represented by each node (binary representation).
is_left    % Indicates for each node if it is a left child.
is_right   % Indicates for each node if it is a right child.
sibling    % Sibling of each node.

end

% Get methods for dependent properties
methods

  function nr_nodes_ = get.nr_nodes(x)
    nr_nodes_ = size(x.children, 1);
  end
   
  function is_leaf_ = get.is_leaf(x)
    is_leaf_ = all(x.children == 0, 2)';
  end

  function lvl_ = get.lvl(x)
    lvl_ = zeros(1, x.nr_nodes);
    ind = htensor.subtree(x.children, 1);
    x_parent = x.parent;
    for ii=ind(2:end)
      lvl_(ii) = lvl_(x_parent(ii)) + 1;
    end
  end
   
  function parent_ = get.parent(x)
    parent_ = zeros(1, x.nr_nodes);
    ind = find(~x.is_leaf);
    parent_(x.children(ind, 1)) = ind;
    parent_(x.children(ind, 2)) = ind;
  end
 
  function dims_ = get.dims(x)
    dims_mat = zeros(1, x.nr_nodes);
    dims_mat(x.dim2ind) = 1:numel(x.dim2ind);
    dims_ = num2cell(dims_mat);
    for ii=x.nr_nodes:-1:1
      if(~x.is_leaf(ii))
        dims_{ii} = [dims_{x.children(ii, 1)}, dims_{x.children(ii, 2)}];
      end
    end
  end

  function bin_dims_ = get.bin_dims(x)
    d = numel(x.dim2ind);
    bin_dims_ = false(x.nr_nodes, d);
    bin_dims_(sub2ind([x.nr_nodes, d], ...
		     x.dim2ind, 1:d)) = true;
    
    for tt=x.nr_nodes:-1:1
      if(~x.is_leaf(tt))
        bin_dims_(tt, :) = bin_dims_(x.children(tt, 1), :) | ...
        bin_dims_(x.children(tt, 2), :);
      end
    end
  end
    
  function is_left_ = get.is_left(x)
    left_nodes = x.children(:, 1);
    left_nodes = left_nodes(left_nodes ~=  0);
    is_left_ = false(x.nr_nodes, 1);
    is_left_(left_nodes) = true;
  end
   
  function is_right_ = get.is_right(x)
    right_nodes = x.children(:, 2);
    right_nodes = right_nodes(right_nodes ~=  0);
    is_right_ = false(x.nr_nodes, 1);
    is_right_(right_nodes) = true;
  end
  
  function sibling_ = get.sibling(x)
    ind = 2:x.nr_nodes;
    child_nodes(ind, :) = x.children(x.parent(ind), :);
    
    x_is_left = x.is_left;
    x_is_right = ~x_is_left; x_is_right(1) = false;
    
    sibling_ = zeros(x.nr_nodes, 1);
    sibling_(x_is_left)  = child_nodes(x_is_left,  2);
    sibling_(x_is_right) = child_nodes(x_is_right, 1);
  end
end


methods( Access = public )

function x = htensor(varargin)
%HTENSOR Construct a tensor in HTD and return htensor object.
%
%   X = HTENSOR(SZ) creates a zero htensor of size SZ with a balanced
%   dimension tree. By default, SZ = [1 1].
%
%   X = HTENSOR(CP) converts a tensor in CP decomposition into an htensor
%   with a balanced dimension tree. CP may be a cell array, or a Tensor
%   Toolbox ktensor object. 
%
%   X = HTENSOR(SZ, TREE_TYPE)
%   X = HTENSOR(CP, TREE_TYPE) allow to specify the type of dimension tree
%   (as described in DEFINE_TREE).
%
%   X = HTENSOR(SZ, CHILDREN, DIM2IND)
%   X = HTENSOR(CP, CHILDREN, DIM2IND) allow to define the dimension tree
%   by the two arrays CHILDREN and DIM2IND. See DEFINE_TREE for a 
%   description of CHILDREN, DIM2IND.
%
%   X = HTENSOR(CHILDREN, DIM2IND, U, B{, IS_ORTHOG}) constructs an htensor
%   directly from the arguments. See DEFINE_TREE for a description
%   of CHILDREN, DIM2IND.
%
%   Examples
%   x = htensor([3 4 5 6]) %<-- creates zero 3x4x5x6 htensor
%   A{1} = rand(3,3); A{2} = rand(4,3); A{3} = rand(5,3); A{4} = rand(6,3);
%   x = htensor(A); %<-- creates a 3x4x5x6 htensor with hierarchical ranks
%                        1,3,...,3 and diagonal transfer tensors
%
%   See also DEFINE_TREE, HTENONES, HTENRANDN, TRUNCATE,
%   (Tensor Toolbox) KTENSOR.


% Default constructor
if(nargin == 0)
  
  x = htensor([1 1]);
  return;
  
elseif(nargin <= 3)
  
  % Conversion constructor from tensor_toolbox class ktensor
  if(isa(varargin{1}, 'ktensor') || isa(varargin{1}, 'cell'))

    % Check for empty cell array
    if(isempty(varargin{1}))
      error(['First argument CP cannot be an empty cell array,' ...
             ' dimension of tensor in CP decomposition must be' ...
             ' at least 1.']);
    end

    % Initialize variables
    if(isa(varargin{1}, 'ktensor'))
      kt = varargin{1};
      d = ndims(kt);
      cpU = kt.U;
      cpU{1} = cpU{1}*diag(kt.lambda);
      rank = length(kt.lambda);
    else
      d = numel(varargin{1});
      cpU = varargin{1};
      rank = size(cpU{1}, 2);
    end

    % Set minimum number of dimensions to 2
    if(d == 1)
      d = 2;
      cpU{2} = ones(1, rank);
    end

    % Check input size and types
    if(~all(cellfun(@(x)(isa(x, 'numeric')), cpU)))
      error('All elements of cell array CP must be numeric.')
    end

    if( any(cellfun('size', cpU, 2) ~= rank) || ...
        any(cellfun('size', cpU, 1) == 0) )
      error(['All arrays in CP must have the same number of columns, ' ...
             'and must contain at least one row.']);
    end

    % Initialize dimension tree
    if(nargin <= 1)
      [x.children, x.dim2ind] = htensor.define_tree(1:d);
    elseif(nargin == 2)
      [x.children, x.dim2ind] = htensor.define_tree(1:d, varargin{2});
    else
      if(isnumeric(varargin{2}) && size(varargin{2}, 1) == 2*d-1)
        x.children = varargin{2};
      else
        error(['Second argument must be either a string option, or' ...
               ' the field children of an htensor']);
      end
      if(isindexvector(varargin{3}) && numel(varargin{3}) == d)
        x.dim2ind = varargin{3};
      else
        error(['Third argument must contain ' ...
               'the field dim2ind of an htensor']);
      end
    end
    
    % Construct cell arrays U and B
    x.U = cell(1, size(x.children, 1));
    x.B = cell(1, size(x.children, 1));
    
    for ii=2:size(x.children, 1)
      if(x.children(ii, 1) == 0)
        x.U{ii} = cpU{x.dim2ind == ii};
      else
        x.B{ii} = diag3d(ones(rank, 1));
      end
    end
    x.B{1} = eye(rank);
    
    x.is_orthog = false;
    
    % Check resulting htensor
    if(nargin == 3)
      check_htensor(x);
    end
    
  elseif( isvector(varargin{1}) )
    % Zero htensor.

    x_size = varargin{1};
    
    % Make sure x_size has at least two entries, creating singleton
    % dimensions if necessary:
    x_size(end+1:2) = 1;
    
    % Check that all elements of x_size are positive integers
    if( any(ceil(x_size) ~= floor(x_size)) || any(x_size <= 0) )
      error('The vector SZ must only contain positive integers.');
    end
    
    d = numel(x_size);
    
    % Initialize dimension tree
    if(nargin==1)
      [x.children, x.dim2ind] = htensor.define_tree(1:length(x_size));
    elseif(nargin == 2)
      [x.children, x.dim2ind] = htensor.define_tree(1:length(x_size), ...
        varargin{2});
    else
      if(isnumeric(varargin{2}) && size(varargin{2}, 1) == 2*d-1)
        x.children = varargin{2};
      else
        error(['Second argument must be either a string option, or' ...
               ' the field children of the new htensor']);
      end
      if(isindexvector(varargin{3}) && numel(varargin{3}) == d)
        x.dim2ind = varargin{3};
      else
        error(['Third argument must contain ' ...
               'the field dim2ind of the new htensor']);
      end
    end
    
    % Fill U and B cell arrays
    x.U = cell(1, size(x.children, 1));
    x.B = cell(1, size(x.children, 1));
    
    for ii=1:size(x.children, 1)
      if(x.children(ii, 1) == 0)
        x.U{ii} = zeros(x_size(x.dim2ind == ii), 1);
        x.B{ii} = [];
      else
        x.U{ii} = [];
        x.B{ii} = 1;
      end
    end
    
    x.is_orthog = false;
    
    % Check resulting htensor
    if(nargin == 3)
      check_htensor(x);
    end
    
  else
    error('Invalid arguments.')
  end

elseif(nargin >= 4)
  % Insert data
  x.children = varargin{1};
  x.dim2ind  = varargin{2};
  x.U        = varargin{3};
  x.B        = varargin{4};
  
  if(nargin >= 5)
    x.is_orthog   = varargin{5};
  else
    x.is_orthog = false;
  end
  
  % Check resulting htensor  
  check_htensor(x);
  
  % Ensure that U and B have correct length
  x.U(end+1:size(x.children, 1)) = {[]};
  x.B(end+1:size(x.children, 1)) = {[]};
  
end

end

% Other public functions.
  AB = apply_mat_to_mat(A, B, p)
  Ax = apply_mat_to_vec(A, x, transp)
  x = cat(ii, x1, x2)
  check_htensor(x)
  x = conj(x)
  ctranspose(x)
  disp(x, name, v)
  disp_all(x, name, v)
  display(x)
  [dofs, dofsB] = dofs(x)
  [z, z_, err, sv] = elem_mult(x, y, opts)
  e = end(x, k, n)
  comp = equal_dimtree(x1, x2)
  names = fieldnames(s)
  y = full(x)
  A = full_mat(A_htd, m, n)
  y = full_block(x, index)
  x = full_leaves(x)
  G = gramians(x)
  G = gramians_nonorthog(x)
  s = innerprod(x, y)
  s = innerprod_mat(x, y, A) 
  x = ipermute(x, order)
  comp = isequal(x1, x2)
  cp = ktensor_approx(x, R, varargin)
  x = minus(x1, x2)
  x = mrdivide(x, a)
  x = mtimes(a, b)
  M = mttkrp(x, U, n)
  d = ndims(x)
  [ndofs, ndofsB] = ndofs(x)
  nrm = norm(x)
  nrm = norm_diff(x, x_full, max_numel)
  x = orthog(x)
  x = permute(x, order)
  plot_sv(x, opts)
  x = plus(x1, x2)
  y = power(x, p)
  r = rank(x, idx)
  sz = size(x, idx)
  sv = singular_values(x, opts)
  x = sparse_leaves(x)
  spy(x, opts)
  y = squeeze(x, dims)
  x = subsasgn(x, s, v)
  out = subsref(x, s)
  transpose(x)
  [x, err, sv] = truncate_std(x, opts)
  [x, err, sv] = truncate_nonorthog(x, opts)
  x_tucker = ttensor(x)
  x = ttm(x, A, varargin)
  x = ttv(x, v, dims)
  z = ttt(x, y, dims_x, dims_y)
  x = uminus(x)
  x = uplus(x)
  x = change_root(x, ind, lr_subtree)
  x = change_dimtree(x, children, dim2ind)
  U = nvecs(x, i, r)
  z = times(x, y)

end

methods( Static, Access = public )

  ind = subtree(children, ii, count)
  [children, dim2ind] = define_tree(dims, opts)
  [U, s] = left_svd_gramian(G)
  [U, s] = left_svd_qr(A)
  [k, err, success] = trunc_rank(s, opts)
  G = gramians_sum(x_cell)
  G = gramians_cp(cp, weights, tree_type)

  [ht, err, sv] = truncate_ltr(x, opts)
  [ht, err, sv] = truncate_rtl(x, opts)
  [x, err, sv] = truncate_sum(x_cell, opts)
  [x, err, sv] = truncate_cp(cp, opts, weights)
     
end

end
