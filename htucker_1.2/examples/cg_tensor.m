function [x, norm_r] = cg_tensor(apply_mat, apply_precond, b, opts, x0)
%Truncated Conjugate Gradient method for htensor.
%
%   [X, NORM_R] = CG_TENSOR(APPLY_MAT, APPLY_PRECOND, B, OPTS) applies the
%   CG method to solve the linear system
%      APPLY_MAT(X) = B,
%   using the preconditioner APPLY_PRECOND. The right-hand side B as well
%   as the solution are htensor objects. The iterates of CG are truncated
%   to low hierarchical rank in each iteration, as specified by OPTS.
%
%   APPLY_MAT(X,OPTS) and APPLY_PRECOND(X,OPTS) are function handles that
%   accept an htensor input X and return a matrix-vector product as
%   htensor.
%
%   Required fields of OPTS:
%   - OPTS.MAXIT       maximal number of iterations of CG method
%   - OPTS.TOL         tolerance for residual norm in CG method
%   - OPTS.MAX_RANK    argument to TRUNCATE
%
%   Optional fields of OPTS:
%   - OPTS.ABS_EPS     argument to TRUNCATE
%   - OPTS.PLOT_CONV   if true, residual norm is plotted during iteration
%   - OPTS.REL_EPS     argument to TRUNCATE
%
%   See also APPLY_LIN_MAT, APPLY_INV_MAT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isfield(opts, 'plot_conv'))
  opts.plot_conv = true;
end

if(nargin == 5)
  x = x0;
else
  x = truncate(0*b, struct('max_rank', 1));
end

r = b - apply_mat(x, opts);
z = apply_precond(r, opts);
p = z;


Ap = apply_mat(p, opts);
pAp = innerprod(p, Ap);

norm_r = zeros(1, opts.maxit);
norm_r(1) = norm(r)/norm(b);

for ii = 2:opts.maxit
  omega = innerprod(r, p)/pAp;
  
  x = x + p*omega;
  x = truncate(x, opts);    
  if(opts.max_rank > 1)
    r = b - apply_mat(x, opts);
  else
    r = b - apply_mat(x);
  end
  r = orthog(r);
  
  norm_r(ii) = norm(r)/norm(b);
  r = truncate(r, opts);
    
  if(opts.plot_conv)
    semilogy(norm_r, 'b');
    hold on;
    drawnow;
  end  
  
  if( norm_r(ii) < opts.tol || ...
      (isfield(opts, 'rel_tol') && ...
       norm_r(ii-1) - norm_r(ii) < opts.rel_tol*norm_r(ii)) )
    break;
  end
  z = apply_precond(r, opts);
  beta = -innerprod(z, Ap)/pAp;
  
  p = z + p*beta;
  p = truncate(p, opts);
    
  Ap = apply_mat(p, opts);
  pAp = innerprod(p, Ap);
end

norm_r = norm_r(1:ii);