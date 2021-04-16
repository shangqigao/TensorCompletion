function x = ipermute(x, order)
%IPERMUTE Inverse permute dimensions of htensor.
%
%   Y = IPERMUTE(X, ORDER) is the inverse of PERMUTE. The dimensions
%   of an htensor X are rearranged such that X = PERMUTE(Y, ORDER). 
%
%   Note that the structure of the dimension tree of X is not modified;
%   only the assigments of the dimensions to the nodes changes.
%
%   Example:
%   x = htenrandn([6,7,8,9]);  %<-- 6x7x8x9 htensor
%   y = permute(x,[3 2 1 4]);  %<-- 8x7x6x9 htensor
%   z = ipermute(y,[3 2 1 4]); %<-- x and z should be identical
%
%   See also IPERMUTE, HTENSOR/PERMUTE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin~=2 )
  error('Two arguments required.');
end
if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end
if(~isindexvector(order))
  error('Second argument ORDER must be a vector of indices.');
end

d = ndims(x);

% Check permutation given by order
if( ~isindexvector(order) || any(order > d) || numel(order) ~= d)
  error(['Invalid permutation, ORDER must be a permutation of 1:' ...
         'ndims(x).']);
end

% Calculate the inverse permutation
inv_order = zeros(1, d);
inv_order(order) = 1:d;

if(any(inv_order == 0))
  error(['Invalid permutation, ORDER must be a permutation of 1:' ...
         'ndims(x).']);
end

% Apply inverse permutation to dims at each node
x.dim2ind = x.dim2ind(inv_order);