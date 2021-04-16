function example_sum_reciproc()
%EXAMPLE_SUM_RECIPROC Apply element-wise reciprocal method to sum of tensors.
%
%   This Matlab function reproduces Example 9.1 from [D. Kressner and
%   C. Tobler. htucker - A Matlab toolbox for the hierarchical Tucker
%   format. Technical report, EPFL, 2012].
%
%   See also HTENSOR, ELEM_RECIPROCAL.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt


close all;

n = 50; d = 4;                

xi = linspace(1, 10, n)';                

xil = xi*ones(1, n^(d-1));
xil = reshape(xil, n*ones(1, d));                

xisum = xil;              
for ii=2:d                
  xisum = xisum + permute(xil, [ii, 2:ii-1, 1, ii+1:d]);              
end
                
x = 1./xisum;

opts.max_rank = 10;  opts.rel_eps = 1e-5;                
x_ht = truncate(x, opts);                
rel_err = norm(x(:) - x_ht(:))/norm(x(:))

plot_sv(x_ht);

M = 25; j = (-M:M);
xMin = d*min(xi);

hst = pi/sqrt(M);
alpha = -2*log(exp(j*hst)+sqrt(1+exp(2*j*hst)))/xMin;
omega = 2*hst./sqrt(1+exp(-2*j*hst))/xMin;

x_cp = cell(1, d);
for ii=1:d
  x_cp{ii} = exp(xi*alpha);
end
x_cp{1} = x_cp{1}*diag(omega);
x_cp = htensor(x_cp);

rel_err = norm(x(:) - x_cp(:))/norm(x(:))

x_trunc_cp = truncate(x_cp, opts);
rel_err = norm(x(:) - x_trunc_cp(:))/norm(x(:))

plot_sv(x_trunc_cp);

% Construction of x using Newton-Schulz iteration

xisum_cp = cell(1, d);
for ii=1:d
  xisum_cp{ii} = ones(n, d);
  xisum_cp{ii}(:, ii) = xi;
end

opts.elem_mult_max_rank = 50;
opts.elem_mult_abs_eps = 1e-2;
opts.max_rank = 50;
opts.rel_eps = 1e-5;

x0 = htenones(size(x)) / ( d*max(xi) );
x_rec = elem_reciprocal(htensor(xisum_cp), opts, x0 );

rel_err = norm(x(:) - x_rec(:))/norm(x(:))

plot_sv(x_rec);

% Quadrature using approximation in HTD

% Construction of quadrature weights 
h = 9/(n-1);
w = 4*ones(n, 1);
w(3:2:end-2) = 2;
w(1) = 1;
w(end) = 1;
w = h/3*w;

w_cell = cell(1, d);

% Inner product between weights and function values by repeated contraction
for ii=1:d
  w_cell{ii} = w;
end

ttv(x_ht, w_cell)