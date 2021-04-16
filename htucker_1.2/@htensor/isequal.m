function comp = isequal(x1, x2)
%ISEQUAL Check whether two htensors are equal.
%
%   ISEQUAL(X1, X2) returns TRUE if all fields of the htensor objects X1
%   and X2 are identical, and FALSE otherwise.
%
%   In particular, htensor objects corresponding to identical tensors 
%   represented with different HTDs are regarded as different objects.
%
%   See also HTENSOR, EQUAL_DIMTREE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires 2 arguments.')
end

if(~isa(x1, 'htensor') || ~isa(x2, 'htensor'))
  error('Both arguments must be of class htensor.');
end

comp = isequal(x1.children, x2.children) & ...
       isequal(x1.dim2ind, x2.dim2ind) & ...
       isequal(x1.U, x2.U) & ...
       isequal(x1.B, x2.B) & ...
       isequal(x1.is_orthog, x2.is_orthog);
