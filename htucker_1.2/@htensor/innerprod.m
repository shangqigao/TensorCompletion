function s = innerprod(x, y)
%INNERPROD Inner product for htensor.
%
%   INNERPROD(X,Y) returns the inner product between two tensors. Both X
%   and Y are assumed to be htensor objects with identical dimension trees
%   and sizes.
%
%   For complex-valued tensors, the complex conjugate of X is used.
%
%   Example:
%   x = htenrandn([3 4 5 6]); y = htenrandn([3 4 5 6]);
%   innerprod(x,y)
%
%   See also HTENSOR, NORM.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin~=2 )
  error('Requires two arguments.');
end
% Check x and y
if(isa(y, 'ktensor'))
  if(~isequal(ndims(x), ndims(y)))
    error('X and Y must be of identical size.');
  end
  y = htensor(y, x.children, x.dim2ind);
elseif(isa(y, 'cell'))
  if(~isequal(ndims(x), numel(y)))
    error('X and Y must be of identical size.');
  end
  y = htensor(y, x.children, x.dim2ind);
end
if(~isa(x, 'htensor') || ~isa(y, 'htensor'))
  error('X and Y must be of class htensor.');
end
if(~equal_dimtree(x, y))
  error('X and Y must have identical dimension trees.');
end
if(~isequal(size(x), size(y)))
  error('X and Y must be of identical size.');
end

M = cell(x.nr_nodes, 1);

x_is_leaf = x.is_leaf;

% Start at leaves, move up the levels
for ii=x.nr_nodes:-1:1
  if(x_is_leaf(ii))
    % M_t = Ux_t' * Uy_t
    M{ii} = full(x.U{ii}'*y.U{ii});
  else
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);
    
    % M_t = Bx_t' * (M_tx kron M_ty) * By_t
    % (interpreting Bx_t, By_t to be in matricized form)
    B_ = ttm(y.B{ii}, { M{ii_left}, M{ii_right} }, [1 2]);
    
    M{ii} = ttt(x.B{ii}, B_, [1 2], [1 2], 3, 3);
    
    % If there is a singleton dimension, M{ii} needs to be
    % reshaped:
    M{ii} = reshape(M{ii}, size(x.B{ii}, 3), size(y.B{ii}, 3));
    
    % Save memory
    M{ii_left} = []; M{ii_right} = [];
  end
end

s = M{1};
