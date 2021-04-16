function [x, err, sv] = truncate_nonorthog(x, opts)
%TRUNCATE_NONORTHOG Truncate htensor to lower-rank htensor.
%
%   Y = TRUNCATE_NONORTHOG(X, OPTS) truncates an htensor X to a lower-rank
%   htensor Y, according to OPTS:
%   OPTS.MAX_RANK (mandatory): maximal hierarchical rank.
%   OPTS.REL_EPS (optional): maximal relative norm-wise error.
%   OPTS.ABS_EPS (optional): maximal absolute norm-wise error.
%
%   Estimates for the overall error and for the error in each node of the
%   dimension tree are displayed if OPTS.DISP_ERRTREE is set to true.
%
%   [Y, ERR, SV] = TRUNCATE_STD(x, opts) also returns the truncation error
%   and singular values in each node of the dimension tree.
%
%   In contrast to TRUNCATE_STD, x is not put in orthogonalized HTD at the
%   start of this method.
%
%   *** This method is only provided for illustration, ***
%   *** TRUNCATE_STD is always faster. However, the    ***
%   *** principle of TRUNCATE_NONORTHOG is used in     ***
%   *** ADD_TRUNCATE and TRUNCATE_CP.                  ***
%
%   See also HTENSOR, TRUNCATE, TRUNCATE_STD, ORTHOG

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires exactly 2 arguments.')
end

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.')
end

if(~isa(opts, 'struct') || ~isfield(opts, 'max_rank') )
  error(['Second argument must be a MATLAB struct with at least one,' ...
	       ' field max_rank.']);
end

% From opts.rel_eps/opts.abs_eps, calculate the error permitted for
% each truncation (there are 2*d - 2 truncations overall).
if(isfield(opts, 'rel_eps'))
  opts.rel_eps = opts.rel_eps/sqrt(2*ndims(x) - 2);
end

if(isfield(opts, 'abs_eps'))
  opts.abs_eps = opts.abs_eps/sqrt(2*ndims(x) - 2);
end

Q_cell = cell(1, x.nr_nodes);
R_cell = cell(1, x.nr_nodes);
for ii=find(x.is_leaf)
  
  % Calculate QR decomposition
  [Q, R] = qr(x.U{ii}, 0);
  
  % Make sure rank doesn't become zero
  if(size(R, 1) == 0)
    Q = ones(size(Q, 1), 1);
    R = ones(1, size(R, 2));
  end
  Q_cell{ii} = Q;
  R_cell{ii} = R;
  
  x.U{ii} = R;
  
end

% Calculate the Gramians of x
G = gramians_nonorthog(x);

% err represents the node-wise truncation errors
err = zeros(x.nr_nodes, 1);

x_is_leaf = x.is_leaf;

sv = cell(1, x.nr_nodes);

% Traverse tree from leaves to root
for ii=x.nr_nodes:-1:2

  if(x_is_leaf(ii))
    
    % Set U{ii} to Q
    x.U{ii} = Q_cell{ii};
    
  else
  
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);
    
    % Multiply Rleft, Rright to B:
    x.B{ii} = ttm(x.B{ii}, {R_cell{ii_left}, R_cell{ii_right}}, [1 2]);
    
    % Matricize B{ii}
    B_mat = matricize(x.B{ii}, [1 2], 3, false);
    
    % Calculate QR-decomposition
    [Q, R] = qr(B_mat, 0);
    
    R_cell{ii} = R;
    
    % Calculate dimensions of "tensor" Q
    tsize_new = size(x.B{ii});
    tsize_new(3) = size(Q, 2);
    
    % Reshape Q to tensor B{ii}
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
  
  % Calculate upper bound and factor c from ||x - x_|| <= c ||x - x_best||
  err_bd = norm(err_); factor = sqrt(2*ndims(x)-3);
  
  fprintf(['\nLower/Upper bound for best approximation error:\n' ...
	   '%e <= ||X - X_best|| <= ||X - X_|| <= %e\n'], ...
	  max(err_bd/factor, max(err)), err_bd);
end
