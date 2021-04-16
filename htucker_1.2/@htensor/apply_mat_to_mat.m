function AB = apply_mat_to_mat(A, B, p)
%APPLY_MAT_TO_MAT Applies an operator in HTD to another operator in HTD.
%
%   AB = APPLY_MAT_TO_MAT(A, B, P) applies the htensor A, which represents
%   a linear operator, to htensor B, which represents another linear
%   operator. The hierarchical ranks of A and B get multiplied.
%
%   A is an htensor of size m_1 p_1 x ... x m_d p_d, while B is an htensor
%   of size p_1 n_1 x ... x p_d n_d. The result AB is an htensor of size
%   m_1 n_1 x ... x m_d n_d. The vector containing the sizes
%   (p_1, ..., p_d) of the inner dimensions is required from the user.
%
%   Examples
%   A = htenrandn([2 4 6 8]); B = htenrandn([4 6 6 2]); p = [2,2,2,2];
%   AB = apply_mat_to_mat(A, B);
%
%   d = 10; n = 100; e = ones(n,1);
%   L = gen_laplace(d,spdiags([e -2*e e], -1:1, n, n));
%   AB = apply_mat_to_mat(L, L, n*ones(1,d));
%
%   See also APPLY_MAT_TO_VEC.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 3)
  error('Requires 3 arguments.');
end

% Check A and B
if(~isa(A, 'htensor'))
  error('First argument must be of class htensor.');
end
if(~isa(B, 'htensor'))
  error('Second argument must be of class htensor.');
end
% Check P
if( ~isindexvector(p) )
  error('p must be a vector of positive integers.');
end

np = size(A);
n = np./p;

mp = size(B);
m = mp./p;

if( any(floor(n) ~= ceil(n)) )
  error('Size mismatch between A and P.');
end

if( any(floor(m) ~= ceil(m)) )
  error('Size mismatch between B and P.');
end

AB = A;

AB_is_leaf = AB.is_leaf;

for ii=1:AB.nr_nodes
  
  if(AB_is_leaf(ii))
    k_A = size(A.U{ii}, 2);
    k_B = size(B.U{ii}, 2);
    mu = find(A.dim2ind == ii);
    
    A_mat = cell(1, k_A);
    for jj=1:k_A
      A_mat{jj} = reshape(A.U{ii}(:, jj), n(mu), p(mu));
    end
    
    B_mat = cell(1, k_B);
    for jj=1:k_B
      B_mat{jj} = reshape(B.U{ii}(:, jj), p(mu), m(mu));
    end
    
    AB_U = cell(k_A, k_B);
    for jj=1:k_A
      for ll=1:k_B
	tmp = A_mat{jj} * B_mat{ll};
	AB_U{ll, jj} = tmp(:);
      end
    end
    AB.U{ii} = cell2mat(transpose(AB_U(:)));
  else
    sz_A = size(A.B{ii});
    sz_A(end+1:3) = 1;
    sz_B = size(B.B{ii});
    sz_B(end+1:3) = 1;
    
    AB.B{ii} = zeros(sz_A.*sz_B);
    for jj=1:sz_A(3)
      for kk=1:sz_B(3)
        AB.B{ii}(:, :, kk+(jj-1)*sz_B(3)) = ...
          kron(A.B{ii}(:, :, jj), B.B{ii}(:, :, kk));
      end
    end
  end
end

AB.is_orthog = false;