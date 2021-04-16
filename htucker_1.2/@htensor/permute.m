function x = permute(x, order)
%PERMUTE Permute dimensions of htensor.
%
%   PERMUTE(X, ORDER) rearranges the dimensions of the htensor X so that
%   they are in the order specified by the vector ORDER.
%
%   Note that the structure of the dimension tree of X is not modified;
%   only the assigments of the dimensions to the nodes changes.
%
%   Example:
%   x = htenrandn([6,7,8,9]); %<-- 6x7x8x9 htensor
%   y = permute(x,[3 2 1 4]); %<-- 8x7x6x9 htensor
%
%   See also PERMUTE, HTENSOR/IPERMUTE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin~=2 )
  error('Two arguments required.');
end
if( ~isa(x, 'htensor') )
  error('First argument must be of class htensor.');
end
if(~isindexvector(order))
  error('Second argument ORDER must be a vector of indices.');
end

d = ndims(x);

% Check permutation given by order
if( ~isindexvector(order) || any(order > d) || length(order) ~= d )
  error(['Invalid permutation, ORDER must be a permutation of 1:' ...
         'ndims(x).']);
end

% Calculate the inverse permutation
inv_order = zeros(1, d);
inv_order(order) = 1:d;

if(any(inv_order == 0))
  error(['Invalid permutation, ORDER must be a permutation of 1:' ...
         'ndims(X).']);
end

% Apply inverse permutation to dim2ind at each node
x.dim2ind = x.dim2ind(order);
