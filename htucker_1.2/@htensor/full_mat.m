function A_full = full_mat(A, m, n)
%FULL_MAT Full matrix represented by an operator in HTD.
%
%   A_FULL = FULL_MAT(A) constructs the matrix represented by htensor A,
%   assuming all matrices are square. Returns an error message if any size
%   of htensor A is not a square number.
%
%   A_FULL = FULL_MAT(A, M, N) constructs the full matrix of size
%   M(1) ... M(D) x N(1) ... N(D) represented by htensor A of size
%   M(1) N(1) x ... x M(D) N(D).
%
%   Example.
%   d = 3; n = 10; e = ones(n,1); L = spdiags([e -2*e e], -1:1, n, n);
%   % 3-dimensional finite-difference Laplace operator
%   A = full_mat(gen_laplace(d, L));

%   See APPLY_MAT_TO_VEC.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin < 1 )
  error('Requires at least one argument.');
end
if( ~isa(A, 'htensor') )
  error('First argument must be of class htensor.');
end

d = ndims(A);
if(nargin == 1)
    n = sqrt(size(A));
    m = n;
elseif(nargin == 2)
    if(~isindexvector(m) || numel(m) ~= d)
        error('M must be index vector with d elements');
    end
    n = m;
elseif(nargin == 3)
    if(~isindexvector(m) || ~isindexvector(n) || ...
       numel(m) ~= d || numel(n) ~= d)
        error('M and N must be index vectors with d elements');
    end
end
if(any(m.*n ~= size(A)))
    error('Sizes given in M and N are not consistent with htensor A');
end

sz = ones(1, 2*d);
sz(1:2:end) = m;
sz(2:2:end) = n;

% Reshape to full tensor of size m1 x n1 x m2 x n2 x ... x md x nd
A_full = reshape(full(A), sz);

% Permute to size m1 x m2 x ... x md x n1 x n2 x ... x nd
A_full = permute(A_full, [1:2:2*d-1, 2:2:2*d]);

% Reshape to matrix
A_full = reshape(A_full, [prod(m), prod(n)]);