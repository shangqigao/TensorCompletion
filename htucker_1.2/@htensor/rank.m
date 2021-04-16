function r = rank(x, ii)
%RANK Hierarchical ranks of htensor.
%
%   RANK(X) returns a (2*NDIMS(X)-1)-vector containing the hierarchical
%   ranks of the htensor X.
%
%   Example:
%   rank(htenones([1 2 3])) %<-- returns [1 1 1 1 1]

%   Internal use:
%   RANK(X,I) returns the ranks for all nodes in the integer array I.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

rk_leaves = cellfun('size', x.U, 2);
rk_nodes = cellfun('size', x.B, 3);
rk_nodes(x.is_leaf) = 0;

r = rk_leaves + rk_nodes;

if nargin > 1
  r = r(ii);
end
