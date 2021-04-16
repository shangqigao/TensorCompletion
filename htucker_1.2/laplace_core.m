function x = laplace_core(d, varargin)
%LAPLACE_CORE Core tensor for Laplace operator.
%
%   X = LAPLACE_CORE(D) generates a 2 x ... x 2 htensor, which is
%   zero everywhere except for
%
%   X(2, 1, ..., 1) = X(1, 2, 1, ..., 1) = ... = X(1, ..., 1, 2) = 1.
%
%   This is the core tensor for the multilinear decomposition of the
%   Laplace operator, and more generally, of any tensor of the form
%  
%   a_1 x b_2 x ... x b_d  +  b_1 x a_2 x ... x b_d  +  ... 
%                     ...  +  b_1 x ... x b_(d-1) x a_d.
%
%   X = LAPLACE_CORE(D, TREE_TYPE)
%   X = LAPLACE_CORE(D, CHILDREN, DIM2IND) allow to choose different types
%   of dimension tree for X, as described in DEFINE_TREE. By default, a
%   balanced dimension tree is chosen.
%
%   Example: % Generates htensor for x(i1, i2, i3, i4) = i1 + i2 + i3 + i4
%   c = laplace_core(4);
%   U = [ones(100, 1), (1:100)'];
%   x = ttm(c, {U, U, U, U});
%
%   See also GEN_LAPLACIAN.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isindexvector(d) || ~isscalar(d) )
  error('First argument must be an integer');
end

% Construct htensor x
x = htensor(2*ones(d, 1), varargin{:});

% Construct U{ii}, B{ii}:
U = cell(1, x.nr_nodes);
B = cell(1, x.nr_nodes);

B_ = dematricize([1 0; 0 1; 0 1; 0 0], [2 2 2], [1 2], 3, false);

x_is_leaf = x.is_leaf;

for ii=2:x.nr_nodes
  
  if( x_is_leaf(ii) )
    U{ii} = speye(2);    
  else
    B{ii} = B_;
  end

end
B{1} = [0 1; 1 0];

x = htensor(x.children, x.dim2ind, U, B, true);