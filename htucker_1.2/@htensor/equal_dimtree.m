function comp = equal_dimtree(x1, x2)
%EQUAL_DIMTREE Compare dimension trees of two htensor objects.
%
%   EQUAL_DIMTREE(X1, X2) returns TRUE if the htensor objects X1 and X2
%   have identical dimension trees, and FALSE otherwise.
%
%   Trees that only become identical after a reordering of the node indices
%   are regarded as different trees.
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires 2 arguments.')
end

if(~isa(x1, 'htensor') || ~isa(x2, 'htensor'))
  error('Both arguments must be of class htensor.');
end

comp = isequal(x1.children, x2.children) & isequal(x1.dim2ind, x2.dim2ind);
