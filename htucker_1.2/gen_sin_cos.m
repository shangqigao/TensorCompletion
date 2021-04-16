function [sin_htd, cos_htd] = gen_sin_cos(t, varargin)
%GEN_SIN_COS Function-valued htensor for sine and cosine.
%
%   [SIN_HTD, COS_HTD] = GEN_SIN_COS(T) generates function-valued tensors
%   for sin(xi_1+...+xi_d) and cos(xi_1+...+xi_d). The cell array T should
%   contain the grid points to be evaluated in each coordinate direction.
%   The returned htensor objects SIN_HTD, COS_HTD have hierarchical ranks 2
%   and satisfy
%   SIN_HTD(i1,..,id) = sin(T{1}(i1) + ... + T{d}(id)),
%   COS_HTD(i1,..,id) = cos(T{1}(i1) + ... + T{d}(id)).
%
%   [SIN_HTD, COS_HTD] = GEN_SIN_COS(T, TREE_TYPE) 
%   [SIN_HTD, COS_HTD] = GEN_SIN_COS(T, CHILDREN, DIM2IND) allow to choose
%   different types of dimension tree for X, as described in DEFINE_TREE.
%   By default, a balanced dimension tree is chosen.
%
%   Example.
%   x = linspace(0,2*pi); t{1} = x; t{2} = x; t{3} = x; t{4} = x;
%   [s, c] = gen_sin_cos(t);
%
%   See also LAPLACE_CORE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if ( nargin < 1 )
  error('Requires at least one argument.');
end

% Check t
if(~iscell(t) || isempty(t) || ~all(cellfun(@isvector, t)) ...
              || ~all(cellfun(@isfloat, t)) )
  error('Argument must be a non-empty cell array containing real/complex vectors');
end

d = numel(t);
% Add zeros if dimension is only one.
if ( d==1 )
  t{2} = 0;
  d = d + 1;
end
sz = cellfun('length', t);

for ii=1:d
  if(size(t{ii}, 1) == 1 && size(t{ii}, 2) > 1)
    t{ii} = t{ii}';
  end
end

% Initialize htensor x
x = htensor(sz, varargin{:});

% Compute U{ii}, B{ii}
U = x.U;
B = x.B;

x_is_leaf = x.is_leaf;

for ii=2:x.nr_nodes
  if(x_is_leaf(ii))
    mu = find(x.dim2ind == ii);
    U{ii} = [sin(t{mu}), cos(t{mu})];
  else
    B{ii} = dematricize([0 -1; 1 0; 1 0; 0 1], [2 2 2], [1 2], 3, false);
  end
end

% Calculate sin_htd
B{1} = [0 1; 1 0];
sin_htd = htensor(x.children, x.dim2ind, U, B);

% Calculate cos_htd
B{1} = [-1 0; 0 1];
cos_htd = htensor(x.children, x.dim2ind, U, B);
