function fHandle = handle_inv_mat(A)
%HANDLE_INV_MAT Function handle to Kronecker structured matrix multiplication.
% 
%   FHANDLE = HANDLE_INV_MAT(A), where A is a matrix or a cell array
%   containing d matrices, returns a function handle FHANDLE. Given an
%   htensor X, a call to FHANDLE(X) performs a matrix-vector
%   multiplication of vec(X) with the Kronecker structured matrix
%
%     inv(A{d}) x ... x inv(A{1}),     if A is a cell array,
%
%     I x ... x I x inv(A),            if A is a matrix.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

  if(isa(A, 'cell'))
    d = length(A);
    for ii=1:d
      A{ii} = @(x)(A{ii}\x);
    end
  else
    d = 1;
    A = @(x)(A\x);
  end

  fHandle = @apply_inv_mat;
  function Ax = apply_inv_mat(x, opts)
    Ax = ttm(x, A, 1:d);
  end
end
