function A = matricize(x, row_dims, col_dims, arg_check)
%MATRICIZE Matricization of (full) tensor.
%
%   A = MATRICIZE(X, ROW_DIMS, COL_DIMS) returns a matricization of X,
%   where ROW_DIMS and COL_DIMS are integer vectors denoting the dimensions
%   of X that are merged into the rows and the columns of the matrix,
%   respectively.
%
%   A = MATRICIZE(X, ROW_DIMS) chooses COL_DIMS automatically to contain
%   all dimensions not contained in ROW_DIMS. The entries of COL_DIMS are
%   arranged in ascending order.
%
%   Note that MATRICIZE(X, 1:ndims(X)) corresponds to X(:).
%
%   Example
%   x = randn(2, 3, 4, 5);
%   A = matricize(x, [2 4]); % returns a 15x8 matrix
%
%   See also DEMATRICIZE.

%   Internal use:
%   A = MATRICIZE(X, ROW_DIMS, COL_DIMS, false) does not check the input
%   arguments.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin == 4 && arg_check == false)
    
  % No argument checking, remove trailing singleton dimensions
  row_dims = row_dims(row_dims <= ndims(x));
  col_dims = col_dims(col_dims <= ndims(x));  
  
else
    
  % Check first two input arguments
  if(nargin <= 1)
    error('Insufficient number of arguments; need at least 2.');
  elseif(~isnumeric(x))
    error('First argument must be numeric.');
  elseif(~isindexvector(row_dims) || ...
    numel(unique(row_dims)) ~= numel(row_dims))
    error('ROW_DIMS must be a positive integer vector without duplicate entries.');
  else
    % Remove trailing singleton dimensions
    row_dims = row_dims(row_dims <= ndims(x));
  end
  
  % Check third (optional) argument
  if(nargin == 2)
    col_dims = setdiff(1:ndims(x), row_dims);
  elseif(isindexvector(col_dims) && numel(unique(col_dims)) == numel(col_dims))
    % Remove trailing singleton dimensions
    col_dims = col_dims(col_dims <= ndims(x));
    
    % Check consistency of input dimensions
    if( ~isequal( sort([row_dims, col_dims]), 1:ndims(x) ) )
      error('ROW_DIMS and COL_DIMS do not match the dimensions of the input tensor.')
    end

  else
    error('COL_DIMS must be a positive integer vector without duplicate entries.')
  end
  
end
  
% Save size of original tensor
sz = size(x);

% Permute dimensions of tensor x
x = permute(x, [row_dims, col_dims]);

% Reshape to matrix
A = reshape(x, prod(sz(row_dims)), prod(sz(col_dims)) );
