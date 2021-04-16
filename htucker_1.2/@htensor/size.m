function sz = size(x, idx)
%SIZE Size of htensor.
%  
%   D = SIZE(X) returns a row vector containing the sizes of the htensor X.
%
%   I = size(X,DIM) returns the size of the dimension specified by the
%   scalar DIM.
%
%   Examples
%   x = htenrandn([3 4 2 1]);
%   size(x)   %<-- returns a length-4 vector
%   size(x,2) %<-- returns 4
%   size(x,5) %<-- ERROR!
%
%   See also HTENSOR, HTENSOR/NDIMS, SIZE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

sz = cellfun('size', x.U(x.dim2ind), 1);

if(nargin > 1)
  if(isindexvector(idx) && isscalar(idx) && idx <= numel(sz))
    sz = sz(idx);
  else
    error(['Second argument must be an integer between 1 and' ...
	         ' %d.'], numel(sz));
  end
end
