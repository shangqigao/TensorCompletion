function x = ttm(x, A, varargin)
%TTM N-mode multiplication of htensor with matrix.
%
%   Y = TTM(X, A, N) computes the N-mode matrix product Y = A o_N X for
%   an htensor X.
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
%   Y = TTM(...,'o') indicates that the orthogonality of the HTD is
%   preserved (only makes sense when all matrices in A are orthogonal). If
%   IS_ORTHOG = true for X then IS_ORTHOG = true for Y.
%
%   Each matrix in A can be replaced by a function handle that performs
%   a matrix-matrix multiplication with the input matrix.
%
%   Examples
%   x = htenrandn([5,3,4,2]);
%   A = rand(4,5); B = rand(4,3); C = rand(3,4); D = rand(3,2); 
%   y = ttm(x, A, 1)             %<-- 1-mode mult. with A
%   y = ttm(x, A, 3, 't')        %<-- 3-mode mult. with A.'
%   y = ttm(x, @fft, 1)          %<-- same as above
%   y = ttm(x, {A,B,C,D})        %<-- 4-way multiplication
%   y = ttm(x, {C,D}, [3 4])  %<-- 3-mode mult. with C, 4-mode mult. with D
%   y = ttm(x, {[],[],C,D})      %<-- same as above
%   y = ttm(x, {A,B,D}, [1 2 4]) %<-- 3-way multiply
%   y = ttm(x, {A,B,C,D}, -3)    %<-- same as above
%   x = orthog(x);               %<-- orthogonalize HTD, x.is_orthog = true
%   Q = orth(rand(4, 4));        %<-- random orthogonal matrix Q
%   y = ttm(x, Q, 3)             %<-- y.is_orthog = false
%   y = ttm(x, Q, 3, 'o')        %<-- y.is_orthog = true
%   y = ttm(x, Q, 3, 't', 'o')   %<-- y.is_orthog = true
%
%   See also TTM, HTENSOR/TTT, HTENSOR/TTV.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt
%
% The structure of this function follows the function ttm from the Tensor
% Toolbox.

% Check number of arguments
if(nargin < 2)
  error('Requires at least 2 arguments.');
end

% Check x
if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

% Check A
% Case of only one matrix: Put it into a cell array.
if (~isa(A, 'cell'))
  if( ~(isfloat(A) && ndims(A) == 2) && ~isa(A, 'function_handle'))
    error(['Second argument must be either a cell array, a' ...
	   ' real/complex matrix or a function handle.']);
  else
    A = {A};
  end
else
  if(~all(cellfun( ...
      @(x)((isfloat(x) && ndims(x) == 2) || ...
	   isa(x, 'function_handle')), A)))
    error(['All elements of cell array A must be either ' ...
	   'matrices or function handles.']);
  end
end

% Default values of optional arguments
transposed = false;
ctransposed = false;
mat_orthog = false;

% Check for additional arguments.
if(nargin >= 3)
  for ii=1:nargin-2
    if(isa(varargin{ii}, 'char') && ~isempty(varargin{ii}))
      if(varargin{ii}(1) == 't')
        transposed = true;
      elseif(varargin{ii}(1) == 'h')
        ctransposed = true;
      elseif(varargin{ii}(1) == 'o')
        mat_orthog = true;
      else
        fprintf('ttm: %dth argument unknown, ignored.\n', ii+2);
      end
    elseif(isindexvector(varargin{ii}) || isindexvector(-varargin{ii}))
      dims = varargin{ii};
    else
      error('Invalid argument DIMS.');
    end
  end
end

% If dims is negative, select all dimensions except those indicated
if(~exist('dims', 'var'))
  
  if(ndims(x) ~= numel(A))
    error(['Invalid argument: number of elements in A must correspond ' ...
           'to the order of tensor X, when DIMS is not given.'])
  end
  dims = 1:ndims(x);
  
elseif( all(dims < 0) )
  
  if(ndims(x) ~= numel(A))
    error(['Invalid argument: number of elements in A must correspond ' ...
           'to the order of tensor X, when negative DIMS are' ...
           ' given.']);
  end
  
  if(any(-dims > ndims(x)))
    error('DIMS contains dimensions outside of the range 1:d.')
  end
  
  dims = setdiff(1:ndims(x), -dims);

  A = A(dims);
else
  if(numel(A) ~= numel(dims))
    error('A must have the same number of elements as DIMS.');
  end
end

% Subsequently compute N-mode matrix products A{1} o_(dims{1}),
% A{2} o_(dims{2}), ...

% Loop over dimensions
for ii=1:length(dims)

  % Ignore empty entries.
  if( isempty(A{ii}) )
    continue;
  end

  % Find node corresponding to dimension dims(ii)
  ind = x.dim2ind(dims(ii));
  
  if(isfloat(A{ii}))
    if(transposed)
      if(size(A{ii}, 1) ~= size(x, dims(ii)))
        error('Matrix dimensions must agree.')
      end
      x.U{ind} = A{ii}.'*x.U{ind};
    elseif(ctransposed)
      if(size(A{ii}, 1) ~= size(x, dims(ii)))
        error('Matrix dimensions must agree.')
      end
      x.U{ind} = A{ii}'*x.U{ind};
    else
      if(size(A{ii}, 2) ~= size(x, dims(ii)))
        error('Matrix dimensions must agree.')
      end
      x.U{ind} = A{ii}*x.U{ind};
    end
  else % A{ii} is a function handle
    x.U{ind} = A{ii}(x.U{ind});
  end 
end

x.is_orthog = x.is_orthog & mat_orthog;
