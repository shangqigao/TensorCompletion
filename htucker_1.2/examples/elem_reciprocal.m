function x =  elem_reciprocal(y, opts, x)
%ELEM_RECIPROCAL Iterative computation of elementwise reciprocal for htensor.
%
%   X = ELEM_RECIPROCAL(Y, OPTS) calculates X = 1./Y iteratively, using the
%   Newton-Schulz iteration:
%      X_(k+1) = X_k + X_k.*(1 - Y.*X_k).
%
%   As each involved element-wise product results in squared hierarchical
%   ranks, truncation is important during this iteration. Truncation may
%   affect the convergence of the Newton-Schulz iteration to a certain
%   extent, depending on the truncation tolerances in OPTS.
%
%   Required fields of OPTS:
%  - MAXIT               Maximal number of Newton-Schultz iterations
%  - MAX_RANK            Maximal truncation rank 
%  - ELEM_MULT_MAX_RANK  Maximal rank during elementwise multiplication
%  - ELEM_MULT_ABS_EPS   Absolute truncation tolerance during elementwise 
%                        multiplication
%
%   X = ELEM_RECIPROCAL(Y, OPTS, X0) uses the htensor X0 as initial guess
%   for the Newton-Schultz iteration. By default, 
%      X_0 = 1/norm(Y)*ALL_ONES,
%   where ALL_ONES is the htensor of all ones.
%   Normally, better choices are given by
%      X_0 = 1/MAX_Y*ALL_ONES,
%   or
%      X_0 = 1/MAX_Y^2*Y,
%   where MAX_Y is a (good) upper bound for the maximal magnitude of all
%   entries in Y.
%
%   See also ELEM_MULT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

opts_elem_mult.max_rank = opts.elem_mult_max_rank;
opts_elem_mult.abs_eps = opts.elem_mult_abs_eps;

if(~isfield(opts, 'verbose'))
  opts.verbose = false;
end

if(~isfield(opts, 'maxit'))
  opts.maxit = 10;
end

all_ones = htenones(size(y), y.children, y.dim2ind);

if(nargin < 3)
  % Choose default initial guess
  x = 1/norm(y)*all_ones;
end

rel_delta = zeros(1, opts.maxit);
for ii=1:opts.maxit
  
  xy = elem_mult(x, y, opts_elem_mult);
  
  delta = truncate(all_ones - xy, opts); 
  
  deltax = elem_mult(delta, x, opts_elem_mult);
  
  x = truncate(x + deltax, opts);
  
  rel_delta(ii) = norm(delta)/norm(all_ones);
  
  if(opts.verbose)
    rel_delta(ii)
    rank_x = rank(x)
    rank_delta = rank(delta)
    
    semilogy(rel_delta);
    ylabel('||y*xk - 1||/||1||')
    drawnow;
  end
  
  if(ii>1 && rel_delta(ii) > rel_delta(ii-1))
    break;
  end
end
