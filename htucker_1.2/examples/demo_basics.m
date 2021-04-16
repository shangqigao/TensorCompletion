function demo_basics()
%DEMO_BASICS Demonstration of basic htensor functionality
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

close all;

disp('> x = htenrandn([3 4 2 3]);');
disp('> y = htenrandn([3 4 2 3]);');
x = htenrandn([3 4 2 3]);
y = htenrandn([3 4 2 3]);

disp('> full(x)');
full(x)
pause;

disp('> rank(x)');
rank(x)
disp('> rank(y)');
rank(y)
pause;

disp('> z = 3*x + y;');
z = 3*x + y;

disp('> rank(z)');
rank(z)
disp('> spy(z)');
spy(z)
pause;

disp('innerprod(z, x)');
innerprod(z, x)

disp('3*innerprod(x, x) + innerprod(y, x)');
3*innerprod(x, x) + innerprod(y, x)
pause;

disp('A = randn(3, 4)');
A = randn(3, 4)

disp('Ax = ttm(x, A, 2);');
Ax = ttm(x, A, 2);
disp('Ax_ = ttm(full(x), A, 2);');
Ax_ = ttm(full(x), A, 2);
disp('norm( Ax(:) - Ax_(:) )');
norm(Ax(:) - Ax_(:))
pause

% use rank = 20, nr_dims = 10 to see that truncate_sum is faster
% than truncate for bigger values.

disp('> d = 5; n = 30; rk = 5;');
disp('> nr_summands = 5; opts.max_rank = 4;');
d = 5; n = 30; rk = 5;
nr_summands = 5; opts.max_rank = 4;

disp('> tensor_sum = htensor(n*ones(1, d));');
disp('> for ii=1:nr_summands');
disp('>   x_cell{ii} = htenrandn(n*ones(1, d), '', rk*ones(2*d-1, 1));');
disp('>   tensor_sum = tensor_sum + x_cell{ii};');
disp('> end');

tensor_sum = htensor(n*ones(1, d));
x_cell = cell(1, nr_summands);
for ii=1:nr_summands
  x_cell{ii} = htenrandn(n*ones(1, d), '', rk*ones(2*d-1, 1));
  tensor_sum = tensor_sum + x_cell{ii};
end

disp(' ');
disp('> spy(tensor_sum);');
disp('> spy(orthog(tensor_sum));');

spy(tensor_sum);
spy(orthog(tensor_sum));
pause;

nrm_tensor_sum = norm(tensor_sum);

disp('> sum_std_trunc = truncate(tensor_sum, opts);');
tic;
sum_std_trunc = truncate(tensor_sum, opts);

fprintf('Elapsed time: %f s\n', toc);
fprintf('Rel. error: %e\n\n', ...
	norm(tensor_sum - sum_std_trunc)/nrm_tensor_sum);


disp('> sum_trunc_sum = htensor.truncate_sum(x_cell, opts);');
tic;
sum_trunc_sum = htensor.truncate_sum(x_cell, opts);

fprintf('Elapsed time: %f s\n', toc);
fprintf('Rel. error: %e\n\n', ...
	norm(tensor_sum - sum_trunc_sum)/nrm_tensor_sum);


disp('> sum_succ_trunc = x_cell{1};');
disp('> for ii=2:nr_summands');
disp('>   sum_succ_trunc = truncate(sum_succ_trunc + x_cell{ii}, opts);');
disp('> end');
tic;
sum_succ_trunc = x_cell{1};
for ii=2:nr_summands
  sum_succ_trunc = truncate(sum_succ_trunc + x_cell{ii}, opts);
end

fprintf('Elapsed time: %f s\n', toc);
fprintf('Rel. error: %e\n\n', ...
	norm(tensor_sum - sum_succ_trunc)/nrm_tensor_sum);


fprintf(['In contrast to truncate, truncate_sum avoids the fill-in\n' ...
      ' shown in the second plot. However, the overhead due to\n' ...
       ' function calls makes it much slower for this small example\n' ...
      ' (see example_truncation.m).\n\n'])

fprintf(['Successive truncation is the fastest for bigger sizes, but\n' ...
      ' may be inaccurate (see example_cancellation.m and\n' ...
      ' example_cancellation2.m)']);
