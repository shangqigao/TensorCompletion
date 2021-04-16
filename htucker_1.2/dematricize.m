function x = dematricize(A, sz, row_dims, col_dims, arg_check)
%DEMATRICIZE Determine (full) tensor from matricization.
%
%   X = DEMATRICIZE(A, SZ, ROW_DIMS, COL_DIMS) returns the tensor X of size
%   SZ. This is the reverse of
%     A = MATRICIZE(X, ROW_DIMS, COL_DIMS)
%   where ROW_DIMS and COL_DIMS contain the dimensions of X that correspond
%   to the rows and to the columns of A, respectively.
%
%   X = DEMATRICIZE(A, SZ, ROW_DIMS) chooses COL_DIMS automatically to
%   contain all dimensions not contained in ROW_DIMS. The entries of
%   COL_DIMS are arranged in ascending order.
%
%   Example:
%   A = randn(15,8);
%   x = dematricize(A, [2 3 4 5], [2 4]) % returns a 2x3x4x5 tensor
%   
%   See also MATRICIZE.

%   Internal use:
%   X = DEMATRICIZE(X, SZ, ROW_DIMS, COL_DIMS, false) does not check the
%   input arguments.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin == 5 && arg_check == false)
  % No argument checking.

  sz(end+1:2) = 1;
  d = numel(sz);
  row_dims = row_dims(row_dims <= d);
  col_dims = col_dims(col_dims <= d);
  
else

  % Check first two input arguments
  if(nargin <= 2)
    error('Insufficient number of arguments; need at least 3.');
  elseif(~isnumeric(A) )
    error('First argument must be numeric.');
  elseif( ~isvector(sz) || any(floor(sz) ~= ceil(sz)) || any(sz < 0) )
    error(['Second argument must be a vector of real non-negative' ...
      ' integers.']);
  elseif(~isindexvector(row_dims) || ...
	 numel(unique(row_dims)) ~= numel(row_dims))
    error('ROW_DIMS must be a positive integer vector without duplicate entries.');
  end
  
  % Insert singleton dimension(s) if numel(sz) < 2
  sz(end+1:2) = 1;
  
  % Remove trailing singleton dimensions
  d = numel(sz);
  row_dims = row_dims(row_dims <= d);
  
  % Check fourth (optional) argument
  if(nargin == 3)
    col_dims = setdiff(1:d, row_dims);
  elseif(isindexvector(col_dims) && numel(unique(col_dims)) == numel(col_dims))
    % Remove trailing singleton dimensions
    col_dims = col_dims(col_dims <= d);

    % Check consistency of input dimensions
    sorted_dims = sort([row_dims, col_dims]);
    if( ~isequal( sorted_dims, 1:d ) )
      error('Dimensions in ROW_DIMS and COL_DIMS do not match.');
    end
  else
    error('COL_DIMS must be a positive integer vector without duplicate entries.')
  end
  
end
  
% Reshape to tensor
x = reshape(A, sz([row_dims, col_dims]) );

% Inverse permute
x = ipermute(x, [row_dims, col_dims]);
