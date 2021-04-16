function y = squeeze(x, dims)
%SQUEEZE Remove singleton dimensions from htensor.
%
%   Y = SQUEEZE(X) for an htensor X returns an htensor Y with the same
%   elements as X, but with all singleton dimensions removed. A singleton
%   dimension is a dimension DIM such that size(X, DIM) == 1.
%
%   Y = SQUEEZE(X, DIMS) only removes the singleton dimensions in the
%   integer array DIMS.
%
%   Note that the returned htensor has at least order 2. In particular, an
%   htensor X of order 2 is unaffected by squeeze.
%
%   Exceptional cases:
%   If all dimensions are singletons (and contained in DIMS if specified).
%
%   Examples
%   squeeze( htenrandn([2 1 3 4]) ) %<-- returns a 2x3x4 htensor
%   squeeze( htenrandn([1 3 1 1]) ) %<-- returns a 3x1 htensor
%   squeeze( htenrandn([1 1 1 1]) ) %<-- returns a scalar
%
%   See also HTENSOR, SQUEEZE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin == 1)
  dims = 1:ndims(x);
elseif(~isindexvector(dims) || max(dims) > ndims(x) || ...
       numel(dims)~=numel(unique(dims)) )
  error('DIMS must be a vector of integers between 1 and ndims(X).')
end

to_squeeze = false(1, ndims(x));
to_squeeze(dims) = true;
if( all(size(x) == 1) && ( length(dims)==ndims(x) ) )
  y = full(x);
  return;
end

% 2-D htensor objects are unaffected by squeeze
if ( ndims(x) == 2 ),
  y = x;
  return;
end

y.size      = size(x);
y.sibling   = x.sibling;
y.parent    = x.parent;
y.dim2ind   = x.dim2ind;
y.is_left   = x.is_left;
y.B         = x.B;
y.U         = x.U;
y.is_leaf   = x.is_leaf;
y.children  = x.children;
y.nr_nodes  = x.nr_nodes;
y.is_orthog = x.is_orthog;

% Remove singleton dimensions
while( length(find(y.size ~= -1)) > 2 )

  % Index of next leaf node to squeeze away
  ind = next_single_node(y, to_squeeze);
  
  % Stop if no singleton dimensions are left
  if(ind == -1)
    break;
  end
  
  % Index of sibling and parent nodes
  ind_sibling   = y.sibling(ind);
  ind_par       = y.parent(ind);
  
 % Dimension corresponding to node ind
  d_ind = find(y.dim2ind == ind);
  
  % Apply U{ii} to parent tensor B, store resulting matrix in tmp
  if(y.is_left(ind))
    tmp = ttm(y.B{ind_par}, y.U{ind}, 1);
    tmp = reshape(tmp, [size(tmp, 2), size(tmp, 3)]);
  else
    tmp = ttm(y.B{ind_par}, y.U{ind}, 2);
    tmp = reshape(tmp, [size(tmp, 1), size(tmp, 3)]);
  end
  
  if(y.is_leaf(ind_sibling))
    % Case sibling node is also a leaf
    
    % Apply sibling node to matrix
    y.U{ind_par} = y.U{ind_sibling}*tmp;
    
    % Set previous B{ind_par} to empty tensor
    y.B{ind_par} = zeros([0 0 0]);
    
    % Update is_leaf
    y.is_leaf(ind_par) = true;
    
    % disconnect ind, ind_sibling from dimension tree
    y.children(ind_par, :)  = [0 0];
    y.parent(ind)  = 0;
    y.parent(ind_sibling) = 0;
    
    % Set dimension d_ind to -1
    y.size(d_ind) = -1;
    
    % Represent d_sibling by ind_par
    y.dim2ind(y.dim2ind == ind_sibling) = ind_par;
    
  else
    % Case sibling node has a subtree
    
    % Apply B{ind_sibling} to B{ind_par}
    y.B{ind_par} = ttm(y.B{ind_sibling}, tmp.', 3);
    
    % Nodes ind and ind_sibling are not used anymore.
    % Change tree to directly connect ind_sibling's
    % children to ind_par.
    y.children(ind_par, :) = y.children(ind_sibling, :);
    y.parent(y.children(ind_sibling, 1)) = ind_par;
    y.parent(y.children(ind_sibling, 2)) = ind_par;   
    
    % Set dimension d_ind to -1
    y.size(d_ind) = -1;
  end
end


% Eliminate unused nodes

% Indices of nodes that are still connected to the tree
new2old_nodes = htensor.subtree(y.children, 1);

% Other direction: new node index for all old nodes
old2new_nodes = zeros(1, y.nr_nodes);
old2new_nodes(new2old_nodes) = 1:length(new2old_nodes);

% Eliminate unused nodes in y.children
y.children = y.children(new2old_nodes, :);
y.dim2ind  = old2new_nodes(y.dim2ind(y.size ~= -1));
y.B        = y.B(new2old_nodes);
y.U        = y.U(new2old_nodes);

% Update node indices in y.children
y.children(y.children ~= 0) = ...
    old2new_nodes(y.children(y.children ~= 0));

y = htensor(y.children, y.dim2ind, y.U, y.B, y.is_orthog);

% In case d == 2 and size == [1 n]: permute
if(ndims(y) == 2 && y.size(1) == 1 && y.size(2) > 1)
  y = permute(y, [2 1]);
end

function ind = next_single_node(x, to_squeeze)
% Returns node indices of a singleton node of x
% If two sibling leaves have singleton dimensions, these nodes are
% returned first.

% All leaf nodes of singleton dimensions
single_nodes = x.dim2ind(x.size == 1 & to_squeeze);

% Parents of two leaf nodes with singleton dimensions
single_node_par = x.parent(single_nodes);
e = ones(numel(single_node_par), 1);
par_two_single_nodes = (sparse(single_node_par, e, e) == 2);
double_nodes = x.children(par_two_single_nodes, 1);

if(~isempty(double_nodes))
  % Return child of parent node, if there is one:
  ind = double_nodes(1);
  
elseif(~isempty(single_nodes))
  % Otherwise, return any singleton node, if one exist
  ind = single_nodes(1);
else
  % No singleton dimensions left
  ind = -1;
end
