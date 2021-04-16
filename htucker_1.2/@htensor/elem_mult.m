function [z, err, sv] = elem_mult(x, y, opts)
%ELEM_MULT Approximate element-by-element multiplication for htensor.
%
%   Z = ELEM_MULT(X, Y, OPTS) performs an approximate element-by-element
%   multiplication of two tensors. Both X and Y are assumed to be htensor
%   objects with identical dimension trees and sizes.
%   In the course of the computation, the resulting htensor Z is truncated
%   to low hierarchical rank according to OPTS:
%   OPTS.MAX_RANK (mandatory): maximal hierarchical rank.
%   OPTS.ABS_EPS (optional): absolute norm-wise error.
%
%   [Z, ERR, SV] = ELEM_MULT(X, Y, OPTS) also returns the truncation error
%   and singular values in each node of the dimension tree.
%
%   Z = ELEM_MULT(X, Y) is identical with Z = X.*Y and performs the
%   element-by-element multiplication exactly, typically at a much higher
%   computational expense.
%   
%   Examples:
%   opts.max_rank = 20; opts.abs_eps = 10^(-10);
%   x = htenrandn([9 8 7 6]); y = htenrandn([9 8 7 6]);
%   z = elem_mult(x, y, opts); 
%
%   x = linspace(0,2*pi); t{1} = x; t{2} = x; t{3} = x; t{4} = x;
%   [s, c] = gen_sin_cos(t);
%   z = elem_mult(s+c, s+c, opts); % Exploits low h. rank of result
%
%   See also TIMES, TRUNCATE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% If opts is not given, calculate exactly
if(nargin == 2)
  z = x.*y;
  return;
end

% Check number of arguments
if(nargin ~= 3)
  error('Requires 2 or 3 arguments.');
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

% Check opts
if(~isa(opts, 'struct') || ~isfield(opts, 'max_rank') )
  error(['Second argument must be a MATLAB struct with at least one,' ...
	       ' field max_rank.']);
end

% opts.rel_eps has no influence:
if(isfield(opts, 'rel_eps'))
  opts = rmfield(opts, 'rel_eps');
end

% Orthogonalize x and y if necessary
x = orthog(x);
y = orthog(y);

% Calculate Gramians
Gx = gramians(x);
Gy = gramians(y);

% Initialize result
z = x;

% err represents the node-wise truncation errors
err = zeros(z.nr_nodes, 1);

U_x = cell(1, z.nr_nodes);
U_y = cell(1, z.nr_nodes);
sv = cell(1, z.nr_nodes);
for ii=2:z.nr_nodes
  
  % Calculate the left singular vectors U_ and singular values s
  % of X_{ii} from the Gramian G{ii}.
  [U_x{ii}, sv_x] = htensor.left_svd_gramian(Gx{ii});
  [U_y{ii}, sv_y] = htensor.left_svd_gramian(Gy{ii});
  
  % Calculate all singular values:
  SV = sv_x*sv_y';
  [sv{ii}, subs] = sort(SV(:), 'descend');
  
  % Calculate rank k to use, and expected error.
  [k, err(ii)] = htensor.trunc_rank(sv{ii}, opts);
  [ind_x, ind_y] = ind2sub(size(SV), subs(1:k));
  
  % Truncate U_x, U_y; S_ = U_x(:, ii) kron U_y(:, ii), ii=1:k
  U_x{ii} = U_x{ii}(:, ind_x);
  U_y{ii} = U_y{ii}(:, ind_y);
  
end

z_is_leaf = z.is_leaf;

for ii=z.nr_nodes:-1:2

  % Apply U_ to node ii:
  if(z_is_leaf(ii))
    Ux = x.U{ii}*U_x{ii};
    Uy = y.U{ii}*U_y{ii};
    z.U{ii} = Ux.*Uy;
  else
    ii_left  = z.children(ii, 1);
    ii_right = z.children(ii, 2);
    
    Bx = ttm(x.B{ii}, {U_x{ii_left}, U_x{ii_right}}, [1 2], 'h');
    Bx = ttm(Bx, U_x{ii}, 3, 't');
    
    By = ttm(y.B{ii}, {U_y{ii_left}, U_y{ii_right}}, [1 2], 'h');
    By = ttm(By, U_y{ii}, 3, 't');
        
    z.B{ii} = Bx.*By;
  end
  
end

ii = 1;
ii_left  = z.children(ii, 1);
ii_right = z.children(ii, 2);

Bx = ttm(x.B{ii}, {U_x{ii_left}, U_x{ii_right}}, [1 2], 'h');
By = ttm(y.B{ii}, {U_y{ii_left}, U_y{ii_right}}, [1 2], 'h');

z.B{ii} = Bx.*By;

z.is_orthog = false;
