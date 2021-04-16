function z = times(x, y)
%TIMES Element-by-element multiplication for htensor.
%
%   Z = TIMES(X, Y) performs an approximate element-by-element
%   multiplication of two tensors. Both X and Y are assumed to be htensor
%   objects with identical dimension trees and sizes.
%
%   Note that this operation is typically quite expensive, as the
%   hierarchical ranks multiply: rank(Z, ii) = rank(X, ii)*rank(Y, ii)
%
%   See also HTENSOR/ELEM_MULT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin~=2 )
  error('Two arguments required.');
end

% Check x and y
if(~isa(x, 'htensor') || ~isa(y, 'htensor'))
  error('X and Y must be of class htensor.');
end
if(~equal_dimtree(x, y))
  error('X and Y must have identical dimension trees.');
end
if(~isequal(size(x), size(y)))
  error('X and Y must be of identical size.');
end

z = x;
z_is_leaf = z.is_leaf;

for ii=1:z.nr_nodes
  
  if(z_is_leaf(ii))
    z.U{ii} = khatrirao_t(x.U{ii}, y.U{ii});
  else
    sz_x = size(x.B{ii});
    sz_x(end+1:3) = 1;
    sz_y = size(y.B{ii});
    sz_y(end+1:3) = 1;

    % Calculate "3-D Kronecker product"
    z.B{ii} = zeros(sz_x.*sz_y);
    for jj=1:size(x.B{ii}, 3)
      for kk=1:size(y.B{ii}, 3)      
        z.B{ii}(:, :, kk+(jj-1)*size(y.B{ii}, 3)) = ...
        kron(x.B{ii}(:, :, jj), y.B{ii}(:, :, kk));
      end
    end
  end
  
end

z.is_orthog = false;