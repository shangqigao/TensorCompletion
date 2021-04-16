function x = reciproc_sum(d, n, a, b, k)
%RECIPROC_SUM Function-valued tensor for 1/(xi_1 + ... + xi_d)
%
%   X = reciproc_sum(D, N, A, B) returns a D-dimensional array X containing
%   the discretization of the function
%            f(xi) = 1 / ( xi_1 + ... + xi_D ),   xi_k in [A,B],
%   on an NxNx...xN uniform grid of the domain [A,B]^D.
%
%   X = reciproc_sum(D, N, A, B, K) returns a cell array containing a CP
%   approximation of tensor rank 2*K+1. This implicitly assumes 0<A<B.
%
%   Examples
%   x = reciproc_sum(3, 100, 0.01, 1);
%   x = reciproc_sum(3, 100, 0.01, 1, 20);

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin < 4 ),
  error('Requires at least 4 input arguments');
end

if(~isindexvector(d) || ~isscalar(d) || ...
   ~isindexvector(n) || ~isscalar(n) || ...
   ~isfloat(a) || ~isscalar(a) || ...
   ~isfloat(b) || ~isscalar(b) || ...
   ( nargin > 4 && ( ~isindexvector(k) || ~isscalar(k) ) ) )
  error(['All arguments must be scalars, and D, N, K must be positive' ...  
         ' integers']);
end

sample1d = linspace(a, b, n).';

if( nargin < 5 ),
  sum_cp = cell(1, d);
  for ii=1:d
    sum_cp{ii} = ones(n, d);
    sum_cp{ii}(:, ii) = sample1d;
  end
  x = 1./(full(htensor(sum_cp)));
else
  % Sinc quadrature rule:
  k_vec = (-k:k);
  hst = pi/sqrt(k);
  t = log(exp(k_vec*hst)+sqrt(1+exp(2*k_vec*hst)));
  w = hst./sqrt(1+exp(-2*k_vec*hst));

  xMin = d*min(sample1d);
  alpha = -2*t/xMin;
  omega = 2*w/xMin;

  % CP decomposition:
  x = cell(1, d);
  for ii=1:d
    x{ii} = exp(sample1d*alpha);
  end

  x{1} = x{1}*diag(omega);
end
