function [x, err, sv] = truncate_cp(cp, opts)
%TRUNCATE_CP Truncate CP tensor to lower-rank htensor.
%
%   Y = TRUNCATE(CP, OPTS) truncates a CP tensor, represented by a cell
%   array or a Tensor Toolbox ktensor object, to a lower-rank htensor Y,
%   according to OPTS:
%   OPTS.MAX_RANK (mandatory): maximal hierarchical rank.
%   OPTS.REL_EPS (optional): maximal relative norm-wise error.
%   OPTS.ABS_EPS (optional): maximal absolute norm-wise error.
%
%   Estimates for the overall error and for the error in each node of the
%   dimension tree are displayed if OPTS.DISP_ERRTREE is set to true.
%
%   OPTS.TREE_TYPE can be used to choose a different type of dimension
%   tree, as described in DEFINE_TREE. By default, a balanced dimension
%   tree is chosen. Alternatively, the dimension tree can be directly
%   defined via OPTS.CHILDREN and OPTS.DIM2IND.
%
%   [Y, ERR, SV] = TRUNCATE_CP(CP, OPTS) also returns the truncation error
%   and singular values in each node of the dimension tree.
%
%   Examples:
%   A{1} = rand(3,5); A{2} = rand(4,5); A{3} = rand(5,5); A{4} = rand(6,5);
%   opts.max_rank = 3;
%   x = htensor.truncate_cp(A,opts); %<-- htucker of rank 3
%
%   opts.max_rank = 20; opts.abs_eps = 0.01; opts.disp_errtree = true;
%   A = reciproc_sum(6, 100, 0.01, 1, 50); %<-- Function-valued CP tensor
%   x = htensor.truncate_cp(A,opts); %<-- htucker with abs. error 0.01
%
%   See also HTENSOR, TRUNCATE, GRAMIANS_CP.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires two arguments.');
end

if(isa(cp, 'ktensor'))
  cp_U = cp.U;
  cp_U{1} = cp_U{1}*diag(cp.lambda);
  d = ndims(cp);
elseif(isa(cp, 'cell'))
  cp_U = cp;
  d = numel(cp);
else
  error('First argument must be a cell array or ktensor object.');
end

if(~all(cellfun(@(x)(isfloat(x)), cp_U)) )
  error('The coefficients of CP must be real or complex matrices.');
end

if(~isa(opts, 'struct') || ~isfield(opts, 'max_rank') )
  error(['Second argument must be a MATLAB struct with at least one,' ...
	       ' field max_rank.']);
end

% From opts.rel_eps/opts.abs_eps, calculate the error permitted for
% each truncation (there are 2*d - 2 truncations overall).
if(isfield(opts, 'rel_eps'))
  opts.rel_eps = opts.rel_eps/sqrt(2*d - 2);
end

if(isfield(opts, 'abs_eps'))
  opts.abs_eps = opts.abs_eps/sqrt(2*d - 2);
end

Q_cell = cell(1, d);
R_cell = cell(1, d);
for ii=1:d
  % Calculate QR decompositions
  [Q, R] = qr(cp_U{ii}, 0);
  
  % Make sure that rank does not become zero
  if(size(R, 1) == 0)
    Q = ones(size(Q, 1), 1);
    R = ones(1, size(R, 2));
  end
  Q_cell{ii} = Q;
  R_cell{ii} = R;
  
  cp_U{ii} = R;  
end

% Calculate the reduced Gramians of x
if(isfield(opts, 'tree_type'))
  G = htensor.gramians_cp(cp_U, opts.tree_type);
  x = htensor(cp_U, opts.tree_type);
elseif(isfield(opts, 'children') && isfield(opts, 'dim2ind'))
  G = htensor.gramians_cp(cp_U, opts.children, opts.dim2ind);
  x = htensor(cp_U, opts.children, opts.dim2ind);
else
  G = htensor.gramians_cp(cp_U);
  x = htensor(cp_U);
end

x.U(x.dim2ind) = Q_cell;

tmp = R_cell;
R_cell = cell(1, x.nr_nodes);
R_cell(x.dim2ind) = tmp;

% err represents the node-wise truncation errors
err = zeros(x.nr_nodes, 1);

x_is_leaf = x.is_leaf;

sv = cell(1, x.nr_nodes);

% Traverse tree from leaves to root
for ii=x.nr_nodes:-1:2
  if(~x_is_leaf(ii))
    
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);
    
    % Multiply Rleft, Rright to B, matricize B{ii}
    B_mat = khatrirao_aux(R_cell{ii_right}, R_cell{ii_left});
    
    % Calculate QR-decomposition
    [Q, R] = qr(B_mat, 0);
    
    R_cell{ii} = R;
    
    % Calculate dimensions of "tensor" Q
    tsize_new = [size(R_cell{ii_left}, 1), ...
		             size(R_cell{ii_right}, 1), size(Q, 2)];
    
    % Reshape Q into tensor B{ii}
    x.B{ii} = dematricize(Q, tsize_new, [1 2], 3, false);
  end
  % Update Gramian
  G{ii} = R_cell{ii}*G{ii}*R_cell{ii}';
  
  % calculate the left singular vectors U_ and singular values s
  % of X_{ii} from the gramian G{ii}.
  [U_, sv{ii}] = htensor.left_svd_gramian(G{ii});
  
  % Calculate rank k to use, and expected error.
  [k, err(ii)] = htensor.trunc_rank(sv{ii}, opts);
  
  % truncate U_
  U_ = U_(:, 1:k);
  
  % Apply U_ to node ii:
  if(x_is_leaf(ii))
    x.U{ii} = x.U{ii}*U_;
  else
    x.B{ii} = ttm(x.B{ii}, U_, 3, 't');
  end
  R_cell{ii} = U_'*R_cell{ii};
end

root_left  = x.children(1, 1);
root_right = x.children(1, 2);

% Multiply Rleft, Rright to B:
x.B{1} = ttm(x.B{1}, {R_cell{root_left}, R_cell{root_right}}, [1 2]);

x.is_orthog = true;

% Display the expected error in each node and overall.
if(isfield(opts, 'disp_errtree') && opts.disp_errtree == true)
  disp(x, 'truncation_error', err);
  
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
  err_ = err; err_(x.children(1, 1)) = 0;

  % Calculate upper bound and c from ||x - x_|| <= c ||x - x_best||
  ind_lvl = find(x.is_leaf);
  factor = sqrt(length(ind_lvl));
  err_bd    = norm(err_(ind_lvl));
  for ii=1:max(x.lvl)
    ind_lvl = find(x.lvl == ii & ~x.is_leaf);
    factor = factor + sqrt(length(ind_lvl));
    err_bd = err_bd + norm(err_(ind_lvl));
  end
  
  fprintf(['\nLower/Upper bound for best approximation error:\n' ...
	   '%e <= ||X - X_best|| <= ||X - X_|| <= %e\n'], ...
	  max(err_bd/factor, max(err)), err_bd);
end
