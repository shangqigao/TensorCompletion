function example_maxest()
%EXAMPLE_MAXEST Example for computing element of maximal absolute value.
%
%   See also HTENSOR, MAXEST.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

close all;

disp('x = htenrandn([40 40 40]);')
x = htenrandn([40 40 40]);

opts.max_rank = 7;
opts.abs_eps = 1e-6;
opts.rel_eps = 1e-6;

opts.elem_mult_max_rank = 30;
opts.elem_mult_abs_eps = 1e-16;

opts.maxit = 30;
% opts.verbose = true;

opts

disp('[max_x, sub] = maxest(x, opts)');
disp('-----------------------------------');
[max_x, sub] = maxest(x, opts);
disp('-----------------------------------');

max_x
sub

[max_exact, ind] = max(abs(x(:)));
[sub1, sub2, sub3] = ind2sub(size(x), ind);

max_exact = x(sub1, sub2, sub3)
sub_exact = [sub1, sub2, sub3]