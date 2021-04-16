function demo_constructor()
%DEMO_CONSTRUCTOR Demonstration of htensor constructors.
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

close all;

disp('Enter return to continue the demo.')
pause

disp('Construct random htensor:')
disp('> x = htenrandn([3 5 6 8])');
x = htenrandn([3 5 6 8])
pause

disp('> norm(x - x)/norm(x)')
nrm = norm(x - x)/norm(x)
disp('> norm(orthog(x - x))/norm(x)')
nrm = norm(orthog(x - x))/norm(x)
disp('> norm(x(:) - x(:))/norm(x(:))')
nrm = norm(x(:) - x(:))/norm(x(:))
pause

if(exist('ktensor.m', 'file'))
  disp('Construct from CP decomposition (using ktensor from Tensor Toolbox):')
  disp('> kt = ktensor(rand(3, 1), {rand(4, 3), rand(5, 3), rand(2, 3), rand(3, 3)})');
  kt = ktensor(rand(3, 1), {rand(4, 3), rand(5, 3), rand(2, 3), rand(3, 3)})
  pause
  disp('> x = htensor(kt)')
  x = htensor(kt)
  pause
  disp('> norm(full(x) - full(kt))')
  nrm = norm(full(x) - full(kt))
  pause
else
  disp('Construct from CP decomposition:')
  disp('> cp_decomp = {rand(4, 3), rand(5, 3), rand(2, 3), rand(3, 3)}');
  cp_decomp = {rand(4, 3), rand(5, 3), rand(2, 3), rand(3, 3)}
  pause
  disp('> x = htensor(cp_decomp)')
  x = htensor(cp_decomp)
  pause
end

disp('Construct by argument:')
disp('> x2 = htensor(x.children, x.dim2ind, x.U, x.B, x.is_orthog);');
x2 = htensor(x.children, x.dim2ind, x.U, x.B, x.is_orthog);
disp('> norm(x(:) - x2(:))/norm(x)')
norm(x(:) - x2(:))/norm(x)
pause

disp('Construct zero tensor:')
disp('> x = htensor([5 2 3 9 5 4])');
x = htensor([5 2 3 9 5 4])
pause

disp('disp(x);');
disp(x);
disp('> x2 = htensor([5 2 3 9 5 4], ''first_separate'');');
x2 = htensor([5 2 3 9 5 4], 'first_separate');
disp('disp(x2);');
disp(x2);
disp('> x3 = htensor([5 2 3 9 5 4], ''TT'');');
x3 = htensor([5 2 3 9 5 4], 'TT');
disp('disp(x3);');
disp(x3);
pause
