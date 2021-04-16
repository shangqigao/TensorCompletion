function [C, U] = gen_invlaplace(A, k, opts)
%GEN_INVLAPLACE htensor for approx. inverse of Laplace-like matrix.
%
%  [C, U] = GEN_INV_LAPLACE(A, k, OPTS) calculates an
%  approximation of the inverse of the Laplace-like matrix
%
%     L = A{1} x I x ... x I + ... + I x ... x I x A{D}.
%
%  The cell array U contains the eigenvectors of A{1}, ..., A{D}, while the
%  htensor C contains an approximation to the eigenvalues of inv(L). This
%  approximation is performed by a combination of sinc-quadrature, using
%  (2k+1) coefficients, and truncation to htensor according to OPTS (see
%  TRUNCATE_CP).
%  
%  C = GEN_INV_LAPLACE(N, k, OPTS), where N is a vector of
%  positive integers, is the special case where A{ii} are the
%  Laplace matrices of size N(ii):
%
%  A = (n+1)^2*spdiags(ones(n, 1)*[-1, 2, -1], [-1 0 1], n, n);
%
%  In this case, U corresponds to the Discrete Sine Transformation,
%  and the values of lambda are explicitly known.
%
% Examples:
%
% x = randn(5, 5, 5); opts.max_rank = 15;
% Lapl = spdiags([-ones(5, 1), 2*ones(5, 1), -ones(5, 1)], ...
% 	      [-1 0 1], 5, 5)*(5+1)^2;
% C = gen_invlaplace([5 5 5], 10, opts);
% Laplx = ttm(x, A, 1) + ttm(x, A, 2) + ttm(x, A, 3);
% x_ = ttm(full(C).*ttm(Laplx, {@dst, @dst, @dst}), ...
%	   {@idst, @idst, @idst});
%
% A = randn(5); A = A*A';
% [C, U] = gen_invlaplace({A, A, A}, 10, opts);
% Laplx = ttm(x, A, 1) + ttm(x, A, 2) + ttm(x, A, 3);
% x_ = ttm(full(C).*ttm(Laplx, U, 'h'), U);

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin ~= 3 ),
  error('Three arguments required.');
end
if( isa(A, 'cell') && ~isempty(A) )

  d = numel(A);  
  n = cellfun('size', A, 1);
  if( ~all(cellfun(@(x)(isfloat(x)), A)) || ...
      ~isequal(n,cellfun('size', A, 2) ) )
    error('A must contain square real or complex matrices.')
  end

  % Diagonalize A{ii}
  U = cell(1, d);
  lambda = cell(1, d);
  for ii=1:d
    [U{ii}, D] = eig(full(A{ii}));
    lambda{ii} = diag(D);
  end
  
elseif(isindexvector(A))
  n = A;
  d = numel(n);
    
  % Initialize eigenvalues of FDM 1d-Laplace matrix
  lambda = cell(1, d);
  for ii=1:d
    lambda{ii} = 4*sin(pi*(1:n(ii))'/(2*(n(ii)+1))).^2*(n(ii)+1)^2;
  end
else
  error(['First argument must be a nonempty cell array, '...
    'or a vector of integers.']);
end

% Quadrature rule:
k_vec = (-k:k);
hst = pi/sqrt(k);
t = log(exp(k_vec*hst)+sqrt(1+exp(2*k_vec*hst)));
w = hst./sqrt(1+exp(-2*k_vec*hst));

lambdaMin = sum(cellfun(@min, lambda));
alpha = -2*t/lambdaMin;
omega = 2*w/lambdaMin;

% CP decomposition:
cell_cols = cell(1, d);
for ii=1:d
  cell_cols{ii} = exp(lambda{ii}*alpha);
end

cell_cols{1} = cell_cols{1}*diag(omega);


% Truncate from CP decomposition to htensor of smaller rank
C = htensor.truncate_cp(cell_cols, opts);
