function demo_elem_reciprocal()
%DEMO_ELEM_RECIPROCAL Demonstration of element-wise reciprocal.
%
%   See also ELEM_RECIPROCAL, HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

close all;

disp('c = laplace_core(4);');
c = laplace_core(4);

disp('U = [ones(100, 1), linspace(1e-3, 1, 100)]''];');
U = [ones(100, 1), linspace(1e-3, 1, 100)'];

disp('x = ttm(c, {U, U, U, U});');
x = ttm(c, {U, U, U, U});

opts.max_rank = 7;
opts.abs_eps = 1e-6;
opts.rel_eps = 1e-6;

opts.elem_mult_max_rank = 30;
opts.elem_mult_abs_eps = 1e-16;

opts.maxit = 30;
opts.verbose = true;

opts

pause;

figure;
disp('inv_x = elem_reciprocal(x, opts);');
disp('-----------------------------------');
inv_x = elem_reciprocal(x, opts);
disp('-----------------------------------');
title('Convergence without knowing maximal entry');

pause;

disp('Upper bound for maximal element of x:');
opts.max_x = 4*max(U(:, 2));

figure;
disp('inv_x = elem_reciprocal(x, opts);');
disp('-----------------------------------');
inv_x = elem_reciprocal(x, opts);
disp('-----------------------------------');
title('Convergence with maximal entry');

disp('plot_sv(inv_x);');
plot_sv(inv_x);
