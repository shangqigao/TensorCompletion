function example_cookies()
%EXAMPLE_COOKIES Apply CG method to a parametric PDE.
%
%   This Matlab function reproduces Example 9.2 from [D. Kressner and
%   C. Tobler. htucker - A Matlab toolbox for the hierarchical Tucker
%   format. Technical report, EPFL, 2012].
%
%   See also HTENSOR, CG_TENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

close all;

load cookies_matrices_2x2

A_handle = handle_lin_mat(A, {[], 0:100, 0:100, 0:100, 0:100});
M_handle = handle_inv_mat(A{1});

e = ones(101, 1);
b_cell = {b, e, e, e, e};
b_tensor = htensor(b_cell);

opts.max_rank = 30;
opts.rel_eps = 1e-10;
opts.maxit = 50;
opts.tol = 0;

[x, norm_r] = cg_tensor(A_handle, M_handle, b_tensor, opts);


x_mean = full(ttv(x, {e,e,e,e}, [2 3 4 5])) / 101^4;

u_mean = U_bd;
u_mean(FreeDofs) = x_mean;

figure;
patch('faces',Mesh.Elements, 'vertices',Mesh.Coordinates, ...
      'FaceVertexCData', u_mean, ...
      'facecolor','interp', 'edgecolor', 'none');
axis tight;

x_diff = x - htensor({x_mean,e,e,e,e});
x_var = diag(full(ttt(x_diff, x_diff, [2 3 4 5]))) / ( 101^4 - 1 );

u_var = U_bd;
u_var(FreeDofs) = x_var;

figure;
patch('faces',Mesh.Elements, 'vertices',Mesh.Coordinates, ...
      'FaceVertexCData', u_var, ...
      'facecolor','interp', 'edgecolor', 'none');

axis tight;