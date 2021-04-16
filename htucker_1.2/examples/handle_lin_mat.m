function fHandle = handle_lin_mat(A, alpha)
%HANDLE_LIN_MAT Function handle to Kronecker structured matrix multiplication.
% 
%   FHANDLE = HANDLE_LIN_MAT(A, ALPHA) for a cell array A containing d
%   matrices and a cell array ALPHA containing d vectors, returns a
%   function handle FHANDLE. Given an htensor X, a call to FHANDLE(X)
%   performs a matrix-vector multiplication of vec(X) with the Kronecker
%   structured matrix
%
%     I x ... x I x A{1} + I x ... x D{2} x A{2} + ... 
%                    ... + D{d} x I x ... x I x A{d},
%
%   where D{i} = diag( ALPHA{i} ) for i = 2, ..., d.
%
%   With an optional second argument OPTS, a call to FHANDLE(X,OPTS)
%   performs truncation after every addition, according to the properties
%   of OPTS as described in TRUNCATE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

  d = length(A);
  m = cellfun('length', alpha);
  
  fHandle = @apply_lin_mat;
  
  function Ax = apply_lin_mat(x, opts)
    Ax = ttm(x, A{1}, 1);
    for jj=2:d
      Ajj_x = ttm(x, A{jj}, 1);
      Ajj_x = ttm(Ajj_x, spdiags(alpha{jj}', 0, m(jj), m(jj)), jj);
      Ax = Ax + Ajj_x;
      
      if(nargin == 2)
         Ax = truncate(Ax, opts);
      end
    end
  end
end

