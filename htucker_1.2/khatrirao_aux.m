function Z = khatrirao_aux(X, Y)
%KHATRIRAO_AUX Khatri-Rao product.
%
%   Z = KHATRIRAO_AUX(X, Y) computes the Khatri-Rao product of two matrices
%   X and Y.
%
%   For an MxN matrix X and an PxN matrix Y, the result is an
%   (M*P)xN matrix Z with
%
%      Z(:, i) = kron(X(:, i), Y(:, i))
%
%   Example:
%   Z = khatrirao_aux(rand(5,2), rand(3,2)) % Returns 15x2 matrix
%
%   See also KHATRIRAO_T, HTENSOR, ELEM_MULT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin < 2)
  error('Insufficient number of input arguments; need at least 2.');
end
if(~isfloat(X) || ~isfloat(Y))
  error('Input arguments must be real or complex arrays.');
end
[m, n] = size(X);
[p, n_] = size(Y);
if(n ~= n_)
  error('Number of columns of X and Y do not match.');
end

Z = zeros(m*p, n);

for i=1:n
  C = Y(:, i)*X(:, i).';
  Z(:, i) = C(:);
end

