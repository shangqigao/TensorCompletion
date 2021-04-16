function Ax = apply_mat_to_vec(A, x, transp)
%APPLY_MAT_TO_VEC Applies an operator in HTD to htensor.
%
%   AX = APPLY_MAT_TO_VEC(A, X) applies the htensor A, which represents an
%   operator, to htensor X. The hierarchcial ranks of A and X get
%   multiplied.
%
%   X is an htensor of size n_1 x ... x n_d, while A is an htensor of size
%   m_1 n_1 x ... x m_d n_d. The result Ax is an htensor of size
%   m_1 x ... x m_d.
%
%   AX = APPLY_MAT_TO_VEC(A, X, 't') applies the transposed of A.
%
%   AX = APPLY_MAT_TO_VEC(A, X, 'h') applies the conjugate transposed of A.
%
%   Example
%   x = htenrandn([1 2 3 4]); A = htenrandn([4 6 6 12]);
%   y = apply_mat_to_vec(A, x);
%
%   See also APPLY_MAT_TO_MAT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin < 2)
  error('Requires at least 2 arguments.');
end

% Check A and x
if(~isa(A, 'htensor'))
  error('First argument must be of class htensor.');
end
if(~isa(x, 'htensor'))
  error('Second argument must be of class htensor.');
end

if( (nargin>=3) && ~ischar(transp) ),
  error('Third argument must be a character.');
end

n = size(x);

nm = size(A);
m = nm./n;

if( any(floor(m) ~= ceil(m)) )
  error('Size mismatch between A and x.')
end

Ax = x;
Ax_is_leaf = Ax.is_leaf;

for ii=1:Ax.nr_nodes
  
  if(Ax_is_leaf(ii))
    k_A = size(A.U{ii}, 2);
    mu = find(x.dim2ind == ii);
    
    AxU = cell(1, k_A);
    for jj=1:k_A
            
      if(nargin == 3 && (transp == 't' || transp == 'h'))
        Ajj = reshape(A.U{ii}(:, jj), n(mu), m(mu));
        if(transp == 't')
          AxU{jj} = Ajj.' * x.U{ii};
        else
          AxU{jj} = Ajj' * x.U{ii};
        end
      else
        Ajj = reshape(A.U{ii}(:, jj), m(mu), n(mu));
        AxU{jj} = Ajj * x.U{ii};
      end
      
    end
    Ax.U{ii} = cell2mat(AxU);
  else
    sz_A = size(A.B{ii});
    sz_A(end+1:3) = 1;
    sz_x = size(x.B{ii});
    sz_x(end+1:3) = 1;

    % "3-D Kronecker product"
    Ax.B{ii} = zeros(sz_A.*sz_x);
    for jj=1:sz_A(3)
      for kk=1:sz_x(3)      
        Ax.B{ii}(:, :, kk+(jj-1)*sz_x(3)) = ...
          kron(A.B{ii}(:, :, jj), x.B{ii}(:, :, kk));
      end
    end
  end
  
end

Ax.is_orthog = false;