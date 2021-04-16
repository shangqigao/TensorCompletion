function [ht, err, sv] = truncate_rtl(x, opts)
%TRUNCATE_RTL Truncate full tensor to htensor, root-to-leaves.
%
%   Y = TRUNCATE_RTL(X, OPTS) truncates a multidimensional array X to an
%   htensor Y, according to OPTS:
%   OPTS.MAX_RANK (mandatory): maximal hierarchical rank.
%   OPTS.REL_EPS (optional): maximal relative norm-wise error.
%   OPTS.ABS_EPS (optional): maximal absolute norm-wise error.
%
%   Alternatively, X may be a Tensor Toolbox tensor.
%
%   *** This method is only provided for illustration, ***
%   *** TRUNCATE_LTR is always faster, and has the     ***
%   *** same error bound.                              ***
%
%   Estimates for the overall error and for the error in each node of the
%   dimension tree are displayed if OPTS.DISP_ERRTREE is set to true.
%
%   OPTS.TREE_TYPE can be used to choose a different type of dimension
%   tree, as described in DEFINE_TREE. By default, a balanced dimension
%   tree is chosen. Alternatively, the dimension tree can be directly
%   defined via OPTS.CHILDREN and OPTS.DIM2IND.
%
%   OPTS.SV controls the computation of singular vectors and values:
%   OPTS.SV = 'svd': From SVD of matricization X_{ii} (default).
%   OPTS.SV = 'gramian': From EIG of Gramian X_{ii}*X_{ii}'.
%
%   [Y, ERR, SV] = TRUNCATE_RTL(X, OPTS) also returns the truncation error
%   and singular values in each node of the dimension tree.
%
%   Examples:
%   x = rand(8,8,8,8,8);
%   opts.max_rank = 3;
%   y = htensor.truncate_rtl(x,opts); %<-- htucker of rank 3
%
%   opts.max_rank = 20; opts.abs_eps = 0.01; opts.disp_errtree = true;
%   A = reciproc_sum(3, 100, 0.01, 1); %<-- Function-valued tensor
%   x = htensor.truncate_rtl(A,opts); %<-- htucker with abs. error 0.01
%
%   See also HTENSOR, TRUNCATE, TRUNCATE_LTR

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires exactly 2 arguments.')
end

if(isa(x, 'tensor'))
  x = double(x);
end

if(~isfloat(x))
  error('First argument must be a real or complex multidimensional array.')
end

if(~isa(opts, 'struct') || ~isfield(opts, 'max_rank') )
  error(['Second argument must be a MATLAB struct with at least one,' ...
	       ' field max_rank.']);
end

% Initialize htensor ht
if(isfield(opts, 'tree_type'))
  ht = htensor(size(x), opts.tree_type);
elseif(isfield(opts, 'children') && isfield(opts, 'dim2ind'))
  ht = htensor(size(x), opts.children, opts.dim2ind);
else
  ht = htensor(size(x));
end

% initialize cells of matrices U and B:
U = cell(1, ht.nr_nodes);
B = cell(1, ht.nr_nodes);

% err represents the node-wise truncation errors
err = zeros(1, ht.nr_nodes);

% From opts.rel_eps/opts.abs_eps, calculate the error permitted for
% each truncation (there are 2*d - 2 truncations overall).
if(isfield(opts, 'rel_eps'))
  opts.rel_eps = opts.rel_eps/sqrt(2*ndims(x) - 2);
end

if(isfield(opts, 'abs_eps'))
  opts.abs_eps = opts.abs_eps/sqrt(2*ndims(x) - 2);
end

ht_is_leaf = ht.is_leaf;
ht_dims = ht.dims;

k = zeros(1, ht.nr_nodes);
sv = cell(1, ht.nr_nodes);

% Traverse tree from leaves to root
for ii=ht.nr_nodes:-1:1
  
  % Matricization of x corresponding to node ii
  x_mat = matricize(x, ht_dims{ii});
  
  % special case of root node: matricization is a vector
  if(ii == 1)
    U{ii} = x_mat;
    k(ii) = 1;
  else
    
    % Calculate left singular vectors U_ and singular values of x_mat
    if(~isfield(opts, 'sv'))
      opts.sv = 'svd';
    end
    
    if(strcmp(opts.sv, 'gramian'))
      [U_, sv{ii}] = htensor.left_svd_gramian(x_mat*x_mat');
    elseif(strcmp(opts.sv, 'svd'))
      [U_, sv{ii}] = htensor.left_svd_qr(x_mat);
    else
      error('Invalid value of OPTS.SV.');
    end
    
    % Calculate rank k to use, and expected error
    [k(ii), err(ii)] = htensor.trunc_rank(sv{ii}, opts);
    
    % Save left singular vectors U for later
    U{ii} = U_(:, 1:k(ii));
  
  end
  
  % Calculate tensor B from U and the Us of the child nodes
  if(~ht_is_leaf(ii))
    
    % Child nodes' indices
    ii_left  = ht.children(ii, 1);
    ii_right = ht.children(ii, 2);
    
    % Number of rows of child node matrices
    n_left = size(U{ii_left}, 1);
    n_right = size(U{ii_right}, 1);
    
    % reshape U{ii} to a k(ii) x n_left x n_right tensor
    U_tensor = dematricize(U{ii}, [n_left, n_right, k(ii)], ...
			   [1 2], 3, false);
    
    % Apply matrices U{ii_left} and U{ii_right} to find
    % transfer tensor B{ii}
    B{ii} = ttm(U_tensor, {U{ii_left}, U{ii_right}}, ...
		     [1 2], 'h');
    
    % Free unused storage
    if(~ht_is_leaf(ii_left))
      U{ii_left} = [];
    end
    if(~ht_is_leaf(ii_right))
      U{ii_right} = [];
    end
  end
  
end

% Free unused storage: 
U{1} = [];

% Call htensor constructor
ht = htensor(ht.children, ht.dim2ind, U, B, false);

% Display the expected error in each node and overall.
if(isfield(opts, 'disp_errtree') && opts.disp_errtree == true)
  disp(ht, 'truncation_error', err);
  
  % We know from theory that
  %
  % ||X - X_best|| <= ||X - X_|| <= err_bd <= factor*||X - X_best||
  %
  % and max(err) <= ||X - X_best||, therefore
  %
  % max(err_bd/factor, max(err)) <= ||X - X_best|| <= ||X - X_|| <= err_bd
  %
  % give upper and lower bounds for the best approximation as well
  % as the truncated version constructed here.
  %
  
  % Count top-level truncation only once
  err_ = err; err_(ht.children(1, 1)) = 0;
  
  % Calculate upper bound and factor c from ||x - x_|| <= c ||x - x_best||
  err_bd = norm(err_); factor = sqrt(2*ndims(ht)-3);
  
  fprintf(['\nLower/Upper bound for best approximation error:\n' ...
	   '%e <= ||X - X_best|| <= ||X - X_|| <= %e\n'], ...
	  max(err_bd/factor, max(err)), err_bd);
end
