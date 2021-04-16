function test_functions()
%TEST_FUNCTIONS Test calls functions of the toolbox
%
% TEST_FUNCTIONS() calls most functions of the toolbox, and is
% intended to detect obvious compatibility problems.
%
% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

matricize(randn([4 6 2 5 3]), [1 4 2]);
matricize(randn([4 6 2 5 3]), [1 4 2], [3 5]);

dematricize(randn(45, 12), [5 3 4 3 3], [1 4 2]);
dematricize(randn(45, 12), [5 3 4 3 3], [1 4 2], [3 5]);

ttm(randn([3 2 4]), randn(5, 2), 2);
ttm(randn([3 2 4]), randn(5, 2), 2);
ttm(randn([3 2 4]), randn(5, 2), 2, 'x');
ttm(randn([3 2 4]), {randn(3), randn(2), randn(4)});
ttm(randn([3 2 4]), {randn(3), randn(2), randn(4)}, -2);
ttm(randn([3 2 4]), {@(x) 2*x, [], @(x) 5*x});

ttt(randn([5 2 4 2]), randn([7 2 3 4]), [4 3], [2 4]);
ttt(randn([2 3 4]), randn([2 4]));
ttt(randn([5 2 4 2]), randn([5 3 4 3]), [1 3]);
ttt(randn([5 2 4 2]), randn([5 2 4 5]), [1 2 3], [4 2 3]);
ttt(randn([5 2 4 5]), randn([5 2 4 1]), [4 2 3], [1 2 3]);


spy3(diag3d(randn(10, 1)));
close;

khatrirao_aux(randn(5, 3), randn(4, 3));

khatrirao_t(randn(3, 4), randn(3, 5));

x = htensor([2 3 5 3]);

y = htensor({randn(3, 2), randn(5, 2), randn(4, 2)}, 'TT');
y = htensor({randn(3, 2), randn(5, 2), randn(4, 2)}, 'first_separate');
y = htensor({randn(3, 2), randn(5, 2), randn(4, 2)}, 'first_pair_separate');
check_htensor(y);

z = htensor({randn(2, 2), randn(3, 2), randn(5, 2), randn(3, 2)}, ...
  x.children, x.dim2ind);

x = htenrandn([2 3 5 3]);
x = htenrandn([2 3 5 3], 'orthog', 5*ones(1, 2*4-1));

y = htenones([2 3 5 3]);
check_htensor(y);

size(x);
ndims(x);
rank(x);
rank(x, 2);

disp(x);
display(x);
disp_all(x);
clc;

ndofs(x);
spy(x);
spy(x, struct('portrait', false, 'title', 'Spider'));
plot_sv(x);
plot_sv(x, struct('portrait', false, 'title', 'SV-Tree'));
x_ = htenrandn([2, 2], 'orthog', [2 2 2]);
x_.B{1} = diag([1, 1e-50]);
plot_sv(x_);
close all;

tmp = 2*x;
tmp = x*4;
tmp = x/3;
conj(x);
tmp = -x;
tmp = +x;

tmp = x(:, 2, 1, end);
tmp = x{1};
tmp = x{4};
tmp = x(1, 2, 1, end);

x.B{1} = ones(size(x.B{1}));
x.U{1} = [];
x.U{4}(1, 1) = 4;
x.B{4} = [];
x{4}(1, 1) = 2;


x = change_root(x, 3);
change_root(x, 1);
change_dimtree(x, x.children(:, [2 1]), x.dim2ind);

x = htenrandn([3 6 2]);
x = change_dimtree(x, x.children(:, [2 1]), x.dim2ind);
x = change_dimtree(x, x.children(:, [2 1]), x.dim2ind);

x = htenrandn([4 6 2 4 6 3 8 2]);
x_ = change_root(x, 5);
change_dimtree(x, x_.children, x_.dim2ind);

x = htenrandn([4 4 4 4 4], 'TT');
x_ = change_dimtree(x, ...
		    [2 3; 4 5; 0 0; 6 7; 0 0; 8 9; 0 0; 0 0; 0 0], ...
		    [8 9 7 5 3]);
change_dimtree(x_, x.children, x.dim2ind);

htensor.subtree(x, 2)

squeeze(htenrandn([3 4 1 2 1]));
squeeze(htenrandn([1 1 1 1 1]));
squeeze(htenrandn([5 1]));
squeeze(htenrandn([1 1 1 5]));

x = htenrandn([4 7 2 3])

y = permute(x, [1 3 2 4]);
x_ = ipermute(y, [1 3 2 4]);
assert(norm(x(:) - x_(:)) == 0);

full(htenrandn([3 2 3 4 2]));
full_block(htenrandn([12, 20, 12 40]), [3 5; 2 4; 9 12; 2 2]);

full_leaves(x);
sparse_leaves(x);

isequal(x, x);

full_mat(htenrandn([4, 9, 9, 16]));
full_mat(htenrandn([4, 9, 9, 16]), [2 3 3 4]);
full_mat(htenrandn([4, 9, 9, 16]), [2 3 1 16], [2 3 9 1]);

orthog(x);
norm(x);
norm(orthog(x));

gramians(x);
singular_values(x);
nvecs(htenrandn([12, 20, 15, 20]), 2, 4);

x = htenrandn([5 3 6 4 4]);
y = htenrandn([5 3 6 4 4]);

tmp = x + y;
tmp = x - y;
innerprod(x,  y);
innerprod(x,  {randn(5, 3), randn(3, 3), randn(6, 3), ...
	       randn(4, 3), randn(4, 3)});
innerprod(change_root(x, 4),  {randn(5, 3), randn(3, 3), randn(6, 3), ...
	       randn(4, 3), randn(4, 3)});

innerprod_mat(htenrandn([2 3 4 5]), htenrandn([2 3 4 5]), ...
	      htenrandn([4 9 16 25]));


norm_diff(x, full(y));
norm_diff(x, full(y), 100);

ttt(x, y, [4 5], [5 4]);
ttt(x, y);
ttt(x, x, 1:5);
ttt(x, y, 3:5);

z = htenrandn([4 5 3]);
w = htenrandn([6 4 4]);
v = htenrandn([8 8 8 8 5 3 4 8]);

ttt(z, x, [2 3 1], [1 2 5]);
ttt(w, x, [1 2 3], 3:5);
ttt(x, w, 3:5, [1 2 3]);
ttt(z, v, [1 2 3], [7 5 6]);

u = htenrandn([8 8 5 3]);
u2 = change_root(u, 2);
ttt(v, u2, [3 4 5 6], [1 2 3 4]);
u3 = change_root(u, 3);
ttt(v, u3, [3 4 5 6], [1 2 3 4]);

ttt(htenrandn(5*ones(9, 1), 'orthog'), ...
    htenrandn(5*ones(3, 1), 'orthog'), ...
    [7 6 5], [3 2 1]);

x = htenrandn([4, 8, 6]);

ttm(x, randn(5, 6), 3);
ttm(x, randn(6, 5), 3, 'h');
ttm(x, randn(6, 5), 3, 't');
ttm(x, orth(randn(10, 6)), 3, 'o');
ttm(x, randn(5, 6), 3, 'x');
ttm(x, {randn(4), randn(8), randn(6)});
ttm(x, {randn(4), randn(8), randn(6)}, -2);
ttm(x, {@(x) 2*x, [], @(x) 5*x});

ttv(x, randn(6, 1), 3);
ttv(x, {[], randn(8, 1), randn(6, 1)});
ttv(x, {randn(4, 1), randn(8, 1), randn(6, 1)}, -3);


x = htenrandn([5 3 6 4 4]);
y = htenrandn([5 3 6 4 4]);

tmp = x.*y;
tmp = x.^2;
elem_mult(x, y, struct('max_rank', 3, ...
		       'rel_eps', 1e-4, ...
		       'abs_eps', 1e-5));
elem_mult(x, y);

cat(2, x, y);

mttkrp(x, {randn(5, 2), randn(3, 2), [], randn(4, 2), randn(4, 2)}, 3);
mttkrp(x, {randn(5, 2), randn(3, 2), [], randn(4, 2), randn(4, 2)}', 3);

A = htenrandn(size(x).^2);
B = htenrandn(size(x).^2);

apply_mat_to_vec(A, x);
apply_mat_to_vec(A, x, 't');
apply_mat_to_vec(A, x, 'h');
apply_mat_to_mat(A, B, size(x));


opts.max_rank = 5;
opts.rel_eps = 1e-8;
opts.disp_errtree = true;
opts.plot_sv = true;

htensor.truncate_rtl(randn(3, 5, 2, 4, 3), opts);

opts.abs_eps = 1e-16;

htensor.truncate_ltr(randn(3, 5, 2, 4, 3), opts);

opts.plot_sv = false;
close all;

truncate_std(x, opts);
truncate_nonorthog(x, opts);

truncate(x, opts);

htensor.truncate_sum({x, y, x}, opts);
htensor.truncate_sum({x, y, x}', opts);
htensor.gramians_sum({x, y, x}');

htensor.truncate_cp({randn(3, 2), randn(5, 2), randn(4, 2)}, opts);

truncate({randn(3, 2), randn(5, 2), randn(4, 2)}, opts);

truncate(randn(3, 2, 4, 5), opts);

opts.tree_type = 'TT';

z = htensor.truncate_rtl(randn(3, 5, 2, 4, 3), opts);
htensor.truncate_ltr(randn(3, 5, 2, 4, 3), opts);
truncate_std(x, opts);
truncate_nonorthog(x, opts);
htensor.truncate_sum({x, y, x}, opts);
htensor.truncate_cp({randn(3, 2), randn(5, 2), randn(4, 2)}, opts);

opts = rmfield(opts, 'tree_type');
opts.children = z.children;
opts.dim2ind = z.dim2ind;
opts.sv = 'gramian';

htensor.truncate_rtl(randn(3, 5, 2, 4, 3), opts);
htensor.truncate_ltr(randn(3, 5, 2, 4, 3), opts);
htensor.truncate_cp({randn(3, 2), randn(5, 2), randn(4, 2), ...
		    randn(7, 2), randn(10, 2)}, opts);


clear all;

laplace_core(10, 'TT');
gen_laplace(10, randn(10), eye(10));
gen_laplace(10, randn(10));
gen_laplace(4, {randn(10), randn(20), randn(30), randn(40)}, ...
	    {randn(10), randn(20), randn(30), randn(40)});

gen_sin_cos({linspace(0, pi), linspace(0, pi), ...
	     linspace(0, pi), linspace(0, pi)});
gen_sin_cos({linspace(0, pi)});

x = reciproc_sum(3, 100, 1, 1e3);
x = reciproc_sum(10, 100, 1, 1e3, 100);

opts.max_rank = 5;
A = randn(10); A = A*A';
gen_invlaplace(repmat({A}, 1, 10), 20, opts);
gen_invlaplace([20 20 20 20], 20, opts);

clc;

% Test the functions that require the Tensor Toolbox, only
% if ktensor is in the MATLAB path.
if(exist('ktensor'))

clear all;

x = htenrandn([5, 5, 5, 5]);

kt = ktensor(rand(3, 1), {randn(2, 3), randn(3, 3), ...
 		    randn(5, 3), randn(3, 3)});

innerprod(htensor(kt), kt);

z = htensor(kt, x.children, x.dim2ind);

opts.max_rank = 2;

htensor.truncate_cp(kt, opts);
truncate(kt, opts);

htensor.truncate_rtl(tensor(randn(3, 5, 2, 4, 3)), opts);
htensor.truncate_ltr(tensor(randn(3, 5, 2, 4, 3)), opts);

truncate(tensor(randn(3, 5, 2, 4, 3)), opts);

kt = ktensor_approx(x, 3, 'maxiters', 5);

tmp = ttensor(x);

clc;

end