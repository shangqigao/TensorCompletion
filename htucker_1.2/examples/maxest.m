function [max_x, ind] = maxest(x, opts)
%MAXEST Approximate element of maximal absolute value.
%
%   [MAX_X, IND] = MAXEST(X, OPTS) aims to compute the element of maximal
%   absolute value for htensor X.
%
%   The computation is performed iteratively by repeated elementwise
%   squaring. As the element-wise product significantly increases
%   hierarchical ranks, truncating is important for this algorithm. On the
%   other hand, too aggressive truncation may harm the reliability of the
%   algorithm. 
%
%   Required fields of OPTS:
%   - ELEM_MULT_MAX_RANK      Option for ELEM_MULT
%   - ELEM_MULT_ABS_EPS       Option for ELEM_MULT
%   - MAX_RANK                Option for TRUNCATE
%   - MAXIT                   Maximal number of iterations
%
%   See also TRUNCATE, ELEM_MULT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

opts_elem_mult.max_rank = opts.elem_mult_max_rank;
opts_elem_mult.abs_eps = opts.elem_mult_abs_eps;

if(~isfield(opts, 'verbose'))
  opts.verbose = false;
end

iter = 1;
y = truncate(x, opts);

% normalize y_trunc
y = orthog(y);
y = y/norm(y);

while( any(rank(y) > 1) && iter < opts.maxit )
  
  [y, err_elmul] = elem_mult(y, y, opts_elem_mult);
  
  % Prevents SVD from returning NaNs in some MATLAB versions:
  for ii=1:x.nr_nodes
    if(x.is_leaf(ii))
      M = y.U{ii};
      M( abs(M) < 1e-20*max(abs(M(:))) ) = 0;
      y.U{ii} = M;
    else
      M = y.B{ii};
      M( abs(M) < 1e-20*max(abs(M(:))) ) = 0;
      y.B{ii} = M;
    end
  end
  
  [y, err_trunc] = truncate(y, opts);
  
  % normalize y
  y = orthog(y);
  y = y/norm(y);
  
  if(opts.verbose)
    rank_y = rank(y)
    
    err_elmul(x.children(1, 1)) = 0;
    err_trunc(x.children(1, 1)) = 0;
    
    err_elmul_trunc = [norm(err_elmul), norm(err_trunc)]
  end
  
  iter = iter + 1;
end

err = zeros(ndims(y), 1);

ind = zeros(1, ndims(y));
for ii=1:ndims(y)
  node_ind = y.dim2ind(ii);
  [val, ind(ii)] = max(abs(y.U{node_ind}));
  err(ii) = abs( abs(val) - 1);
end

max_x = full_block(x, [ind', ind']);
