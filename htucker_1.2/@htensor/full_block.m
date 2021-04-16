function y = full_block(x, index)
%FULL_BLOCK Return subblock of htensor as a (full) tensor.
%
%   Y = FULL_BLOCK(X, INDEX) converts a subblock of htensor X into a (full)
%   multidimensional array Y. The Dx2-matrix INDEX gives the start and end
%   index in each dimension.
%
%   Example:
%   x = htenrandn([2 4 3 5]);
%   y = full_block(x, [1 2; 1 1; 2 3; 1 5])
%   z = full(x(1:2, 1, 2:3, :));   % should be equal to y
%
%   See also HTENSOR, FULL, SUBSREF.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Check x
if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

% Check index
if(~isnumeric(index) || any(size(index) ~= [ndims(x), 2]) || ...
   any(ceil(index(:)) ~= floor(index(:))) || any(index(:) <= 0) )
  error(['Second argument INDEX must be a Dx2-matrix of indices,' ...
         ' where D = ndims(X).']);
end

if ( any( index(:,2)'>size(x) ) || any( index(:,1)>index(:,2) )  ),
  error('Indices out of range.');
end

% Go through all leaves
for ii=find(x.is_leaf)
  % Restrict leaf matrix to the elements given by index
  mu = find(x.dim2ind == ii);
  x.U{ii} = x.U{ii}(index(mu, 1):index(mu, 2), :);
end

% convert x to multdimensional array y
y = full(x);
