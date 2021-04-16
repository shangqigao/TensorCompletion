function x = htenrandn(sz, orthog, k, varargin)
%HTENRANDN Random htensor.
%
% The entries of U_t, B_t are normally distributed pseudo-random numbers.
%
%   X = HTENRANDN(SZ) returns an htensor of size SZ with
%     - hierarchical ranks randomly distributed between 3 and 6;
%     - all elements of the leaf bases U_t and transfer tensors B_t
%       normally distributed pseudo-random numbers;
%     - a balanced dimension tree.
%
%   X = HTENRANDN(SZ, ORTHOG) returns an htensor in orthogonalized HTD if
%   ORTHOG='orthog'.
%
%   X = HTENRANDN(SZ, ORTHOG, K) allows to prescribe the hierarchical ranks
%   for each node in a (2*ndims(X)-1)-vector K. In case of an infeasible
%   choice, the actual hierarchical ranks of X will be lower.
%
%   X = HTENRANDN(SZ, ORTHOG, K, TREE_TYPE) allows to specify the type of
%   dimension tree (as described in DEFINE_TREE).
%
%   X = HTENRANDN(SZ, ORTHOG, K, CHILDREN, DIM2IND) allows to define the
%   dimension tree by the two arrays CHILDREN and DIM2IND.
%
%   Examples
%   htenrandn([3 4 5 6])           % <-- returns 3x4x5x6 htensor
%   htenrandn([3 4 5 6],'orthog')  % <-- returns orthog. 3x4x5x6 htensor
%   X = htenrandn([3 4 5 6],'',[2 2 2 2 2 2 2]) % <-- all ranks 2
%   rank(X) % <-- returns [1 2 2 2 2 2 2]; the root must have rank 1

%   See also HTENSOR, HTENSOR/DEFINE_TREE, HTENONES, RANDN,
%   (Tensor Toolbox TENRAND, SPTENRAND).

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin < 1)
  error('Requires at least 1 argument.');
elseif(nargin == 1)
  orthog = '';
end

% Check sz
if( ~isindexvector(sz) )
  error('SZ must be a vector of positive integers.');
end
if( isempty(sz) )
  error('SZ must not be empty.');
end

% Check orthog
if(~ischar(orthog))
  error('Second argument ORTHOG must be a string.');
end

% Construct htensor x with zero entries
x = htensor(sz, varargin{:});
sz = size(x);

% Random k if none is assigned
if(nargin <= 2 || isempty(k))
  k = floor(4*rand(x.nr_nodes, 1))+3;
end

% Check k
if( ~isindexvector(k) || length(k) ~= x.nr_nodes)
  error('K must be a vector of 2*d-1 positive integers.');
end

% Make sure that the rank in dimension ii is not greater than the
% size in dimension ii.
x_dim2ind = x.dim2ind;

for ii=1:ndims(x)
  if( k(x_dim2ind(ii)) > sz(ii) )
    k(x_dim2ind(ii)) = sz(ii);
  end
end

% Set root rank to 1
k(1) = 1;

x_is_leaf = x.is_leaf;

U = cell(1, x.nr_nodes);
B = cell(1, x.nr_nodes);

for ii=1:x.nr_nodes

  if(x_is_leaf(ii))
    
    if(strcmp(orthog, 'orthog'))
      % Initialize orthogonal U{ii}
      U{ii} = orth(randn(size(x.U{ii}, 1), k(ii)));
    else
      % Initialize non-orthogonal U{ii}
      U{ii} = randn(size(x.U{ii}, 1), k(ii));
    end
    
  else
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);
    
    % Construct random tensor B{ii}
    B{ii} = randn([k(ii_left), k(ii_right), k(ii)]);
    
    if(strcmp(orthog, 'orthog'))
      % Matricize B{ii}
      B_mat = matricize(B{ii}, [1 2], 3, false);
      
      % Calculate orthonormal basis
      Q_mat = orth(B_mat);
      
      % Reshape Q_mat to 3d-array B{ii}
      B{ii} = dematricize(Q_mat, size(B{ii}), [1 2], 3, false);
    end
    
  end
end

x = htensor(x.children, x.dim2ind, U, B, strcmp(orthog, 'orthog'));
