function [U, s] = left_svd_qr(A)
%LEFT_SVD_QR Left singular vectors and values.
%
%   [U,S] = LEFT_SVD_QR(A) computes the left singular vectors U and
%   singular values of the matrix A.
%
%   See also LEFT_SVD_GRAMIAN.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isfloat(A) || ndims(A) ~= 2)
  error('First argument must be a real or complex matrix.');
end

if(size(A, 1) > size(A, 2))
  [Q, R] = qr(A, 0);
  [U, S, V] = svd(R, 0);
  U = Q*U;
else
  R = qr(A');
  R = triu(R);
  R = R(1:size(R, 2), :);
  [U, S, V] = svd(R', 0);
end
s = diag(S);
