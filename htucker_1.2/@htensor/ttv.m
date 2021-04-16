function x = ttv(x, v, dims)
%TTV Tensor-times-vector for htensor.
%
%   Y = TTV(X, V, N) for a column vector V computes the N-mode product of
%   htensor X with V. Note that ndims(Y) = ndims(X) - 1 because the N-th
%   dimension is a singleton and removed.
%
%   Y = TTV(X, V) for a cell array V, containing ndims(X) column vectors,
%   subsequently computes all N-mode products with the vectors contained in
%   V. For an empty entry in V, nothing is done in the corresponding
%   dimension.
%
%   Y = TTV(X, V, DIMS) for a cell array V computes the sequence of N-mode
%   htensor-vector products along the modes specified by DIMS.
%
%   Y = TTV(X, V, EXCLUDE_DIMS) for a cell array V and a vector
%   EXCLUDE_DIMS of negative integers computes all htensor-vector products
%   except at dimensions i for which -i is contained in EXCLUDE_DIMS.
%
%   Examples
%   x = htenrandn([5,3,4,2]);
%   a = rand(5,1); b = rand(3,1); c = rand(4,1); d = rand(2,1);
%   y = ttv(x, a, 1)           %<-- mode-1 multiplication of x with a
%   y = ttv(X, {a,b,c,d})      %<-- multiplication in all modes
%   y = ttv(x, {c,d}, [3 4])   %<-- mode-3 and mode-4 multiplication
%   Y = ttv(X, {a,b,d}, -3)    %<-- same as above
%   Y = ttv(X, {a,b,c,d}, -3)  %<-- same as above
%
%   See also HTENSOR, HTENSOR/TTT, HTENSOR/TTM, HTENSOR/SQUEEZE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt
%
% The structure of this function follows the function ttv from the Tensor
% Toolbox.

% Check number of arguments
if(nargin < 2)
  error('Requires at least 2 arguments.');
end

% Check class of x
if(~isa(x, 'htensor'))
  error('First argument X must be of class htensor.');
end

% Default values of optional arguments
if(nargin == 2)
  dims = 1:ndims(x);
end

% Case of only one vector: Put it into a cell array.
if (isa(v, 'cell') && ~isempty(v) )
  if ( ~all( cellfun(@(x)(isfloat(x)),v ) )  ),
    error('All elements of cell array V must be real or complex vectors.');
  end
elseif(isvector(v) && isfloat(v) )
  v = {v};
else
  error('V must be a vector or a nonempty cell array of real or complex vectors.')
end

% Check dims:
if(~isindexvector(dims) && ~isindexvector(-dims))
  error('DIMS must be a vector of indices.');
end

% If dims is negative, select all dimensions except those indicated
if( all(dims < 0) )
  dims = setdiff(1:ndims(x), -dims);
end

% Check that dims are valid
if( any(dims > ndims(x)) || numel(dims) ~= numel(unique(dims)) )
  error('Invalid argument DIMS.');
end

%  Loop over cell array v
for ii=1:length(dims)
  % Ignore empty entries.
  if( isempty(v{ii}) )
    continue;
  end

  ind = x.dim2ind(dims(ii));
  if(size(v{ii}, 1) ~= size(x, dims(ii)))
    error('Tensor size must match vector lengths.')
  end
  x.U{ind} = v{ii}.'*x.U{ind};
end

x.is_orthog = false;
x = squeeze(x, dims);
