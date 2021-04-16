function x = minus(x1, x2)
%MINUS Binary subtraction for htensor.
%
%   X = MINUS(X1,X2) computes X1-X2 for two htensor objects X1 and X2 of the
%   same size.
%
%   See also HTENSOR, HTENSOR/PLUS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin~=2 )
  error('Two arguments required.');
end
% Check type of x1, x2
if( ~isa(x1, 'htensor') || ~isa(x2, 'htensor') )
  error('Both arguments must be of class htensor.');
end
% Check compatibility of dimension trees.
if(~equal_dimtree(x1, x2))
  error('Dimension trees of X1 and X2 differ.');
end
% Check sizes
if(~isequal(size(x1), size(x2)))
  error('X1 and X2 must be of identical size.')
end

x = plus(x1, -x2);