function [ht, err, sv] = truncate_ltr(x, opts)
%TRUNCATE_LTR Truncate full tensor to htensor, leaves-to-root.
%
%   Y = TRUNCATE_LTR(X, OPTS) truncates a multidimensional array X to an
%   htensor Y, according to OPTS:
%   OPTS.MAX_RANK (mandatory): maximal hierarchical rank.
%   OPTS.REL_EPS (optional): maximal relative norm-wise error.
%   OPTS.ABS_EPS (optional): maximal absolute norm-wise error.
%
%   Alternatively, X may be a Tensor Toolbox tensor.
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
%   [Y, ERR, SV] = TRUNCATE_LTR(X, OPTS) also returns the truncation error
%   and singular values in each node of the dimension tree.
%
%   Examples:
%   x = rand(8,8,8,8,8);
%   opts.max_rank = 3;
%   y = htensor.truncate_ltr(x,opts); %<-- htucker of rank 3
%
%   opts.max_rank = 20; opts.abs_eps = 0.01; opts.disp_errtree = true;
%   A = reciproc_sum(3, 100, 0.01, 1); %<-- Function-valued tensor
%   x = htensor.truncate_ltr(A,opts); %<-- htucker with abs. error 0.01
%
%   See also HTENSOR, TRUNCATE, TRUNCATE_RTL.

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

% Initialize htensor t
if(isfield(opts, 'tree_type'))
  ht = htensor(size(x), opts.tree_type);
elseif(isfield(opts, 'children') && isfield(opts, 'dim2ind'))
  ht = htensor(size(x), opts.children, opts.dim2ind);
else
  ht = htensor(size(x));
end

if(~isfield(opts, 'sv'))
  opts.sv = 'svd';
end

% initialize cells of matrices U and B:
U = cell(1, ht.nr_nodes);
B = cell(1, ht.nr_nodes);

% Make a temporary copy of tensor x
x_ = full(x);

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

% Traverse leaves of the dimension tree
sv = cell(1, ht.nr_nodes);
k = zeros(1, ht.nr_nodes);
for ii=find(ht_is_leaf)
  % Matricization of x corresponding to node ii
  x_mat = matricize(x, ht_dims{ii}, ...
                    [1:ht_dims{ii}-1, ht_dims{ii}+1:ndims(x)], false);
  
  % Calculate left singular vectors U_ and singular values of x_mat
  if(strcmp(opts.sv, 'gramian'))
    [U_, sv{ii}] = htensor.left_svd_gramian(x_mat*x_mat');
  elseif(strcmp(opts.sv, 'svd'))
    [U_, sv{ii}] = htensor.left_svd_qr(x_mat);
  else
    error('Invalid value of OPTS.SV.');
  end
  
  % Calculate rank k to use, and expected error.
  [k(ii), err(ii)] = htensor.trunc_rank(sv{ii}, opts);
  
  % Save left singular vectors U for later
  U{ii} = U_(:, 1:k(ii));
  
  % Reduce tensor x_ along this dimension
  x_ = ttm(x_, U{ii}, ht_dims{ii}, 'h');
end

% Set x to be the reduced tensor x_
x = x_;

ht_lvl = ht.lvl;

% Go through all levels from leaves to root node
for lvl_iter = max(ht_lvl):-1:0
  % Go through all nodes at given level
  for ii=find(ht_lvl == lvl_iter)
    
    % Leaves have already been treated
    if(ht_is_leaf(ii))
      continue;
    end
    
    % Matricization of x corresponding to node ii
    x_mat = matricize(x, ht_dims{ii});
    
    % special case root node: matricization is a vector
    if(ii == 1)   
      U_ = x_mat;
      k(ii) = 1;
    else
      
      % Calculate left singular vectors U_ and singular values of x_mat
      if(strcmp(opts.sv, 'gramian'))
        [U_, sv{ii}] = htensor.left_svd_gramian(x_mat*x_mat');
      elseif(strcmp(opts.sv, 'svd'))
        [U_, sv{ii}] = htensor.left_svd_qr(x_mat);
      else
        error('Invalid argument OPTS.SV.');
      end
      
      % Calculate rank k to use, and expected error.
      [k(ii), err(ii)] = htensor.trunc_rank(sv{ii}, opts);
      
      % Cut U_ after first k columns
      U_ = U_(:, 1:k(ii));
    end
    
    % Child nodes' indices
    ii_left  = ht.children(ii, 1);
    ii_right = ht.children(ii, 2);
    
    % reshape B{ii} from matrix U_ to a 
    % k(ii) x k(ii_left) x k(ii_right) tensor, 
    B{ii} = dematricize(U_, [k(ii_left), k(ii_right), k(ii)], ...
                        [1 2], 3, false);
			  
    % Reduce tensor x_ along dimensions x.dims{ii}; this will
    % change the number of dimensions of x_:
    
    % Matricization of x_, making dims{ii} the row dimensions
    x_mat_ = matricize(x_, ht_dims{ii});
    
    % calculate B{ii}'*x_mat_
    U_x_mat = U_'*x_mat_;
    
    % Instead of eliminating one of the dimensions, just set
    % it to be a singleton, to keep the dimension order consistent
    tsize_red = size(x_); tsize_red(end+1:ndims(ht)) = 1;
    tsize_red(ht_dims{ii_left }(1)) = k(ii);
    tsize_red(ht_dims{ii_right}(1)) = 1;
    
    % Reshape x_mat_ to tensor x_
    x_ = dematricize(U_x_mat, tsize_red, ht_dims{ii});
  end
  
  % Set x to be the reduced tensor x_  
  x = x_;
end

% Call htensor constructor
ht = htensor(ht.children, ht.dim2ind, U, B, true);

% Display the estimated errors
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
  
  % Count top-level truncation only once
  err_ = err; err_(ht.children(1, 1)) = 0;

  % Calculate upper bound and factor c from ||x - x_|| <= c ||x - x_best||
  err_bd = norm(err_); factor = sqrt(2*ndims(ht)-3);
  
  fprintf(['\nLower/Upper bound for best approximation error:\n' ...
	   '%e <= ||X - X_best|| <= ||X - X_|| <= %e\n'], ...
	  max(err_bd/factor, max(err)), err_bd);
end
