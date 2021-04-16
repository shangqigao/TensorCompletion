function L = gen_laplace(d, A, M, varargin)
%GEN_LAPLACE htensor for Laplace-like matrix.
%
%   L = GEN_LAPLACE(d, A, M) returns an htensor representing the matrix
%       M{D} x    ...            x M{2} x A{1} 
%     + M{D} x    ...     x M{3} x A{2} x M{1} 
%     +                   ... 
%     + A{D} x M{D-1} x   ...    x M{2} x M{1},
%   where x denotes the Kronecker product.
%
%   If A is a single matrix instead of a cell array, then identical
%   matrices are used in each dimension: A{1} = ... = A{d} = A. An
%   analogous statement holds for M. In addition, if no third argument is
%   given, all matrices M are assumed to identity matrices.
%
%   The operator L is represented as an htensor of size
%      n_1^2 x ... x n_d^2,
%   by vectorising all matrices, and storing these vectors as (sparse) leaf
%   matrices.
%
%   L = GEN_LAPLACE(d, A, M, TREE_TYPE)
%   L = GEN_LAPLACE(d, A, M, CHILDREN, DIM2IND)
%   allow to choose different types of dimension tree for X, as described
%   in DEFINE_TREE. By default, a balanced dimension tree is chosen.
%
%   Example.
%   d = 10; n = 100; e = ones(n,1); A = spdiags([e -2*e e], -1:1, n, n);
%   % 10-dimensional finite-difference Laplace operator
%   L = gen_laplace(d, A);
%
%   See also LAPLACE_CORE, GEN_INVLAPLACE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin < 2)
  error('Requires at least 2 arguments.');
end
if( ~isindexvector(d) || ~isscalar(d) || d < 2 )
  error('First argument must be an integer greater than one.');
end

% Check cell array A
if(iscell(A))
  if length(A) < d,
    error('Cell array A must have at least D elements.');
  end
  n = cellfun('size', A, 1);
  if( ~all(cellfun(@(x)(isfloat(x)), A)) || ...
      ~all(cellfun(@(x)(ndims(x)==2), A)) || ...
      ~isequal(n,cellfun('size', A, 2) ) )
    error('All coefficients of A must be real or complex square matrices.');
  end
elseif(isfloat(A) && ndims(A) == 2)
  if ( size(A,1) ~= size(A,2) )
    error('A must be a square matrix.');
  end
  n = size(A, 1)*ones(1, d);
  A = {A};
  for ii=2:d
    A{ii} = A{1};
  end
else
  error('Second argument must be a cell array or a real/complex matrix.');
end

% Check or initialize cell array M
if(nargin == 2 || isempty(M))
  M = cell(1, d);
  for ii=1:d
    M{ii} = speye(n(ii));
  end
elseif(iscell(M))
  if length(M) < d,
    error('Cell array M must have at least D elements.');
  end
  if( ~all(cellfun(@(x)(isfloat(x)), M)) || ...
      ~all(cellfun(@(x)(ndims(x)==2), M)) || ...
      ~isequal(n,cellfun('size', M, 1) ) || ...
      ~isequal(n,cellfun('size', M, 2) ) )
    error(['All coefficients of M must be real or complex matrices,', ...
          'matching the sizes of the coefficients in A.']);
  end
elseif(isfloat(M) && ndims(M) == 2)
  m = size(M, 1)*ones(1, d);
  if( ( size(M,1)~=size(M,2) ) || ~isequal(m,n) )
    error('M must be a square matrix matching the size of A.');
  end
  M = {M};
  for ii=2:d
    M{ii} = M{1};
  end
else
  error('Third argument must be a cell array or real/complex matrix.');
end

% Initialize htensor Lapl
Lapl = htensor(1:d, varargin{:});

U = cell(1, Lapl.nr_nodes);
B = cell(1, Lapl.nr_nodes);

% Core tensor (see LAPLACE_CORE.m)
B_ = dematricize([1 0; 0 1; 0 1; 0 0], [2 2 2], [1 2], 3, false);

Lapl_is_leaf = Lapl.is_leaf;

% Compute U{ii}, B{ii}
for ii=2:Lapl.nr_nodes
  if( Lapl_is_leaf(ii) )
    mu = find(Lapl.dim2ind == ii);
    U{ii} = [M{mu}(:), A{mu}(:)];
  else
    B{ii} = B_;
  end
end

B{1} = [0 1; 1 0];

% Construct htensor L
L = htensor(Lapl.children, Lapl.dim2ind, U, B, false);