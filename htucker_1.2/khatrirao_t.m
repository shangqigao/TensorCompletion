function Z = khatrirao_t(X, Y)
%KHATRIRAO_T Transposed Khatri-Rao product.
%
%   Z = KHATRIRAO_T(X, Y) computes the transposed Khatri-Rao product of two
%   matrices X and Y.
%
%   For an NxM matrix X and an NxP matrix Y, the result is an Nx(M*P)
%   matrix Z with
%
%      Z(i, :) = kron(X(i, :), Y(i, :))
%
%   Example:
%   Z = khatrirao_t(rand(2,5), rand(2,3)) % Returns 2x15 matrix
%
%   See also KHATRIRAO_AUX, HTENSOR, ELEM_MULT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin < 2)
  error('Insufficient number of input arguments; need at least 2.');
end
if(~isfloat(X) || ~isfloat(Y))
  error('Input arguments must be real or complex arrays.');
end
[n, m] = size(X);
[n_, p] = size(Y);
if(n ~= n_)
  error('Number of rows of X and Y do not match.');
end

Z = zeros(n, m*p);

for i=1:m
  for j=1:p
    Z(:, j + p*(i-1)) = X(:, i).*Y(:, j);
  end
end