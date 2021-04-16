function d = ndims(x)
%NDIMS Order (number of dimensions) of htensor.
%
%   NDIMS(X) returns the order of htensor X.
%
%   Examples
%   X = htenrandn([4,3]);   ndims(X) %<-- Returns 2
%   X = htenrandn([4 3 1]); ndims(X) %<-- Returns 3
%
%   See also HTENSOR, NDIMS, HTENSOR/SIZE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

d = numel(x.dim2ind);