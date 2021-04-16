function [ndofs, ndofsB] = ndofs(x)
%NDOFS Number of degrees of freedom in htensor.
%
%   [NDOFS, NDOFSB] = NDOFS(X) returns the number of degrees of freedom in
%   all nodes (NDOFS) and only in non-leaf nodes (NDOFSB).
%
%   See also NUMEL.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

ndofsB = sum(cellfun(@numel, x.B));
ndofsU = sum(cellfun(@numel, x.U));

ndofs = ndofsU + ndofsB;