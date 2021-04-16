function x = ttm(x, A, varargin)
%TTM N-mode multiplication of (full) tensor with matrix
%
%   Y = TTM(X, A, N) computes the N-mode matrix product Y = A o_N X for
%   a multi-dimensional array X.
%   The arguments must satisfy size(X, N) == size(A, 2), and the resulting
%   array Y satisfies size(Y, N) == size(A, 1).
%
%   Y = TTM(X, A) for a cell array A subsequently computes the N-mode
%   products with the matrices contained in A:
%     Y = ... A{2} o_DIMS(2) A{1} o_DIMS(1) X
%   For an empty entry in A, nothing is done in the corresponding dimension.
%
%   Y = TTM(X, A, DIMS) for a cell array A and a vector DIMS of positive
%   integers computes the N-mode matrix products along the dimensions
%   specified in DIMS:
%     Y = ... A{2} o_DIMS(2) A{1} o_DIMS(1) X
%   Each dimension must appear at most once in DIMS.
%
%   Y = TTM(X, A, EXCLUDE_DIMS) for a cell array A and a vector
%   EXCLUDE_DIMS of negative integers computes all N-mode matrix products
%   except at dimensions i for which -i is contained in EXCLUDE_DIMS.
%
%   Y = TTM(...,'t') performs the same multiplications as above, but with
%   the (non-conjugate) transposed matrices.
%
%   Y = TTM(...,'h') performs the same multiplications as above, but with
%   the conjugate transposed matrices.
%
%   Each matrix in A can be replaced by a function handle that performs
%   a matrix-matrix multiplication with the input matrix.
%
%   Examples
%   x = rand([5,3,4,2]);
%   A = rand(4,5); B = rand(4,3); C = rand(3,4); D = rand(3,2); 
%   y = ttm(x, A, 1)          %<-- 1-mode mult. with A
%   y = ttm(x, A, 3, 't')     %<-- 3-mode mult. with A.'
%   y = ttm(x, @fft, 1)       %<-- same as above
%   y = ttm(x, {A,B,C,D})     %<-- 4-way multiplication
%   y = ttm(x, {C,D}, [3 4])  %<-- 3-mode mult. with C, 4-mode mult. with D
%   y = ttm(x, {[],[],C,D})   %<-- same as above
%   y = ttm(x, {A,B,D}, [1 2 4])  %<-- 3-way multiply
%   y = ttm(x, {A,B,C,D}, -3) %<-- same as above
%
%   See also TTT, HTENSOR/TTM, HTENSOR/TTT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt
%
% The structure of this function follows the function ttm from the Tensor
% Toolbox.

% Check number of arguments
if(nargin < 2)
  error('Requires at least two arguments.');
end

if(~isfloat(x))
  error('First argument must be a MATLAB multidimensional real or complex array.');
end

% Case of only one matrix: Put it into a cell array.
if (~isa(A, 'cell'))
  if(~isfloat(A) && ~isa(A, 'function_handle'))
    error(['Second argument must be either a cell array, a' ...
           ' real/complex matrix, or a function handle.']);
  else
    A = {A};
  end
else
  if(~all(cellfun( ...
      @(x)(isfloat(x) || isa(x, 'function_handle')), A)))
    error(['All elements of cell array A must be either ' ...
           'real/complex matrices or function handles.']);
  end
end

% Default values of optional arguments
transposed = false;
ctransposed = false;

% Check for additional arguments.
for ii=1:numel(varargin)
  if(isa(varargin{ii}, 'char') && ~isempty(varargin{ii}))
    if(varargin{ii}(1) == 't')
      transposed = true;
    elseif(varargin{ii}(1) == 'h')
      ctransposed = true;
    else
      fprintf('%dth argument in ttm unknown, ignored.\n', ii+2);
    end
  elseif(isindexvector(varargin{ii}) || isindexvector(-varargin{ii}))
    dims = varargin{ii};
  else
    error('Invalid argument in ttm.');
  end
end

if(~exist('dims', 'var'))
  % Set dims to 1, 2, 3, ...
  d = max(ndims(x), numel(A));
  dims = 1:numel(A);
elseif(all(dims < 0))
  % Set dims to all dimensions not excluded.
  d = max(ndims(x), numel(A));
  if(any(-dims > d))
    error('DIMS contains entries out of range.')
  end
  dims = setdiff(1:numel(A), -dims);
  A = A(dims);
else
  d = max(ndims(x), max(dims));
  if(numel(A) ~= numel(dims))
    error('A must have the same number of elements as DIMS.');
  end
end

% Subsequently compute N-mode matrix products A{1} o_(dims{1}),
% A{2} o_(dims{2}), ...

% Calculate size (add singleton dimensions if A contains more entries)
sz = size(x);
sz(end+1:d) = 1;

% Loop over dimensions
for ii=1:length(dims)

  % Ignore empty entries.
  if( isempty(A{ii}) )
    continue;
  end
  
  % Matricize tensor x.
  compl_dims = [1:dims(ii)-1, dims(ii)+1:d];
  
  if(dims(ii) == 1 || isa(A{ii}, 'function_handle'))
    X = matricize(x, dims(ii), compl_dims, false);
    transX = false;
  else
    X = matricize(x, compl_dims, dims(ii), false);
    transX = true;
  end
  
  % Apply matrix A{ii}
  if(isfloat(A{ii}))
    if(transposed)
      if(size(A{ii}, 1) ~= sz(dims(ii)))
        error('Matrix dimensions must agree.')
      end
      if(~transX)
        X = A{ii}.'*X;
      else
        X = X*A{ii};
      end
    elseif(ctransposed)
      if(size(A{ii}, 1) ~= sz(dims(ii)))
        error('Matrix dimensions must agree.')
      end
      if(~transX)
        X = A{ii}'*X;
      else
        X = X*conj(A{ii});
      end
    else
      if(size(A{ii}, 2) ~= sz(dims(ii)))
        error('Matrix dimensions must agree.')
      end
      if(~transX)
        X = A{ii}*X;
      else
        X = X*A{ii}.';
      end
    end
  elseif(isa(A{ii}, 'function_handle'))
    X = A{ii}(X);
  else
    error(['A{ii} is of incompatible type; should be a matrix or' ...
      ' function handle.'])
  end
  
  % Update size of X
  if(transX == false)
    sz(dims(ii)) = size(X, 1);
  else
    sz(dims(ii)) = size(X, 2);
  end

  % Dematricize X
  if(transX == false)
    x = dematricize(X, sz, dims(ii), compl_dims, false);
  else
    x = dematricize(X, sz, compl_dims, dims(ii), false);
  end
end
