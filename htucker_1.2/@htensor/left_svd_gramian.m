function [U, s] = left_svd_gramian(G)
%LEFT_SVD_GRAMIAN Left singular vectors and values from Gramian.
%
%   LEFT_SVD_QR(G) computes the left singular vectors and singular values
%   of a matrix A from its Gramian G = A*A'. It is assumed that G is
%   symmetric or Hermitian.
%
%   See also LEFT_SVD_QR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isfloat(G) || ndims(G)~=2 || size(G, 1) ~= size(G, 2))
  error('G must be a square real or complex matrix.');
end

% Symmetrize matrix G
G = (G + G')/2;

% Compute symmetric eigenvalue decomposition
[U, S] = svd(G);

% Compute singular values from eigenvalues of the gramian
s = sqrt(abs(diag(S)));

% Sort singular values and vectors
[s, idx] = sort(s, 'descend');
U = U(:, idx);
