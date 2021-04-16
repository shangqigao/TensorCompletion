function [x, err, sv] = truncate_std(x, opts)
%TRUNCATE_STD Truncate htensor to lower-rank htensor.
%
%   Y = TRUNCATE_STD(X, OPTS) truncates an htensor X to a lower-rank
%   htensor Y, according to OPTS:
%   OPTS.MAX_RANK (mandatory): maximal hierarchical rank.
%   OPTS.REL_EPS (optional): maximal relative norm-wise error.
%   OPTS.ABS_EPS (optional): maximal absolute norm-wise error.
%
%   Estimates for the overall error and for the error in each node of the
%   dimension tree are displayed if OPTS.DISP_ERRTREE is set to true.
%
%   [Y, ERR, SV] = TRUNCATE_STD(X, OPTS) also returns the truncation error
%   and singular values in each node of the dimension tree.
%
%   Example:
%   x = htenrandn([3 4 5 6],'',[1 10 10 10 10 10 10]);
%   opts.max_rank = 3;
%   y = truncate_std(x,opts)
%
%   See also HTENSOR, TRUNCATE, ORTHOG, GRAMIANS

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

% Orthogonalize x (does nothing if x is already orthogonal)
x = orthog(x);

% Calculate the gramians of x
G = gramians(x);

% err represents the node-wise truncation errors
err = zeros(x.nr_nodes, 1);

x_is_leaf = x.is_leaf;
x_is_left = x.is_left;
x_parent  = x.parent;

sv = cell(1, x.nr_nodes);

% Go from leaves to root (though the order does not matter)
for ii=x.nr_nodes:-1:2

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
  
  % Parent node
  ii_par = x_parent(ii);
  
  % Apply U_ to parent node
  if(x_is_left(ii))
    x.B{ii_par} = ttm(x.B{ii_par}, U_, 1, 'h');
  else
    x.B{ii_par} = ttm(x.B{ii_par}, U_, 2, 'h');
  end
  
end

x.is_orthog = false;

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
