function s = innerprod_mat(x, y, A)
%INNERPROD_MAT Weighted inner product for htensor.
%
%   S = INNERPROD_MAT(X, Y, A) returns the inner product between two
%   htensor objects X and Y with respect to a matrix represented by
%   htensor. All htensor objects A, X, and Y are assumed to have
%   identical dimension trees and matching sizes.
%
%   For complex-valued tensors, the complex conjugate of X is used.
%
%   This is a more efficient way of calculating 
%      innerprod(X, apply_mat_to_vec(A, Y));
%
%   Example:
%   x = htenrandn([2 3 4 5]); y = htenrandn([2 3 4 5]);
%   A = htenrandn([4 9 16 25]);
%   innerprod_mat(x, y, A);
%
%   See also APPLY_MAT_TO_VEC, INNERPROD.

if( nargin~=3 )
  error('Requires three arguments.');
end

% Check x, y, and A.
if(~isa(x, 'htensor') || ~isa(y, 'htensor') || ~isa(A, 'htensor'))
  error('X, Y, and A must be of class htensor.');
end
if(~equal_dimtree(x, y) || ~equal_dimtree(x, A))
  error('X, Y, and A must have identical dimension trees.');
end
if(~isequal(size(x), size(y)) || ~isequal(size(x).^2, size(A)))
  error('Sizes of X, Y, and A do not match.')
end

M = cell(1, x.nr_nodes);
n = size(x); 
x_is_leaf = x.is_leaf;

for jj=x.nr_nodes:-1:2
  
  if(x_is_leaf(jj))
    k1 = size(x.U{jj}, 2);
    k2 = size(y.U{jj}, 2);
    s = size(A.U{jj}, 2);
    mu = find(x.dim2ind == jj);
    
    M{jj} = zeros(k1, s, k2);
    
    for ii=1:s
      Aii = reshape(A.U{jj}(:, ii), n(mu), n(mu));
      M{jj}(:, ii, :) = x.U{jj}' * Aii * y.U{jj};
    end
    
  else
        
    jj1 = x.children(jj, 1);
    jj2 = x.children(jj, 2);
    
    [k1l, k1r, k1p] = size(x.B{jj});
    [k2l, k2r, k2p] = size(y.B{jj});
    [sl , sr , sp ] = size( A.B{jj});    
    
    B_1 = ttt(M{jj1}, x.B{jj}, 1, 1, [2 3], [2 3]);
    B_1 = reshape(B_1, [sl, k2l, k1r, k1p]);
    B_2 = ttt(M{jj2}, y.B{jj}, 3, 2, [1 2], [1 3]);
    B_2 = reshape(B_2, [k1r, sr, k2l, k2p]);
    
    B_ = ttt(B_1, B_2, [2 3], [3 1], [1 4], [2 4]);
    B_ = reshape(B_, [sl, k1p, sr, k2p]);
    M_ = ttt(B_, A.B{jj}, [1 3], [1 2], [2 4], 3);
    M_ = reshape(M_, [k1p, k2p, sp]);
    M{jj} = permute(M_, [1 3 2]);
    M{jj} = reshape(M{jj}, [k1p, sp, k2p]);
    
    M{jj1} = [];
    M{jj2} = [];
    
  end
  
end

jj = 1;
jj1 = x.children(jj, 1);
jj2 = x.children(jj, 2);

M_ = ttm(M{jj2}, {x.B{1}, A.B{1}, y.B{1}});
s = M{jj1}(:)'*M_(:);
