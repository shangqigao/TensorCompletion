function M = mttkrp(x, U, n)
%MTTKRP Building block for approximating htensor by ktensor.
%
%   V = MTTKRP(X, U, N) returns the matrix product of the N-mode
%   matricization of X with the Khatri-Rao product of all entries (except
%   the Nth entry) in the cell array U.
%
%   See also KTENSOR_APPROX, (Tensor Toolbox) CP_ALS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 3)
  error('Requires exactly 3 arguments.');
end

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

if(~isa(U, 'cell') || ~all(cellfun( ...
    @(x)((isfloat(x) && ndims(x)==2) || isa(x, 'function_handle')), U)))
  error(['Second argument U must be a cell array of real/complex matrices or' ...
	 ' function handles.']);
end

if(size(U, 2) == 1)
  U = U.';
end

if(~(isindexvector(n) && isscalar(n)))
  error('Third argument N must be a positive integer.')
end

% Order and dimensions of x
d = ndims(x);
sz = size(x);

% Check order of input
if(ndims(x) ~= length(U))
  error(['Number of elements in cell array U must be equal to the' ...
	 ' order of tensor X.']);
end

m1 = cellfun('size', U, 1); m1_red = m1([1:n-1,n+1:d]);
m2 = cellfun('size', U, 2); m2_red = m2([1:n-1,n+1:d]);

% Check number of rows of U
if( any(sz([1:n-1,n+1:d]) ~= m1_red) )
  error('Dimension mismatch.');
end

% Check number of columns of U
if( any(diff(m2_red ~= 0) ) )
  error(['All matrices in cell array U must have the same number' ...
	 ' of columns.'])
end

% Apply U{ii}' to x in all dimensions except n
x = ttm(x, U, -n, 't');

% The following is a slow version of what is calculated in the
% subsequent code:
% M = zeros(sz(n), m);
% for jj=1:m
%   ind = jj*ones(d, 2);
%   ind(n, :) = [1, sz(n)];
%   M(:, jj) = full_block(x, ind);
% end

% Change node corresponding to mode n, s.t. it becomes the
% left child of the root node
x = change_root(x, x.dim2ind(n), 'left');

x_is_leaf = x.is_leaf;

% Loop through all nodes except root node and its left child:
for ii=x.nr_nodes:-1:3
  if(~x_is_leaf(ii))
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);

    x.U{ii} = khatrirao_t(x.U{ii_right}, x.U{ii_left}) * ...
	      matricize(x.B{ii}, [1 2], 3, false);
  end
end

ii = 1;
ii_left  = x.children(ii, 1);
ii_right = x.children(ii, 2);

M = x.U{ii_left}*x.B{ii}'*x.U{ii_right}.';
