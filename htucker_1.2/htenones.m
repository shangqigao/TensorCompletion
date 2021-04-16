function x = htenones(sz, varargin)
%HTENONES htensor with all elements one
%
%   X = HTENONES(SZ) returns an htensor of size SZ with all elements one
%   and a balanced dimension tree.
%
%   X = HTENONES(SZ, TREE_TYPE) allows to specify the type of dimension
%   tree (as described in DEFINE_TREE).
%
%   X = HTENONES(SZ, CHILDREN, DIM2IND) allows to define the dimension tree
%   by the two arrays CHILDREN and DIM2IND.
%
%   Example
%   htenones([2 3 4]) %<-- 2x3x4 tensor
%
%   See also HTENSOR, HTENSOR/DEFINE_TREE, HTENRANDN, ONES.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin < 1)
  error('Requires at least 1 argument.');
end

% Check sz
if( ~isindexvector(sz) )
  error('SZ must be a vector of positive integers.');
end
if( isempty(sz) )
  error('SZ must not be empty.');
end

% Construct htensor x with zero entries
x = htensor(sz, varargin{:});

sz = size(x);
x_is_leaf = x.is_leaf;

% Construct htensor x with zero entries
U = cell(1, numel(sz));
B = cell(1, numel(sz));
for ii=1:x.nr_nodes
  if(x_is_leaf(ii))
    % Initialize U{ii}
    U{ii} = ones(sz(x.dim2ind == ii), 1);
  else
    % Set B{ii}
    B{ii} = 1;
  end
end

x = htensor(x.children, x.dim2ind, U, B);
