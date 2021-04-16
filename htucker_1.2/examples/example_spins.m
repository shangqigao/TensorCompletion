function example_spins(d)
%EXAMPLE_SPINS Demonstration of operator-HTD for 1D spin system.
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

close all;

fprintf(['Consider a 1D system of spins with nearest neighbor interaction\n' ...
         'and periodic boundary conditions (example from Huckle et\n' ...
         'al., Computations in Quantum Tensor Networks):\n'])
fprintf(['  P x P x I x ... x I + I x P x P x I x ... x I + ...\n' ...
         '  + I x ... x P x P + P x I x ... x I x P.\n\n']); 

if(nargin == 0)
  d = 8;
end

fprintf('d = %d\n\n', d);

disp('Construct core tensor only:');
disp('for ii=1:d');
disp('  core_kt{ii} = [ones(1, d); zeros(1, d)];');
disp('  core_kt{ii}(:, ii) = [0; 1];');
disp('  core_kt{ii}(:, mod(ii, d)+1) = [0; 1];');
disp('end');
disp('core_htd = htensor(core_kt);'); disp(' ');
core_kt = cell(1, d);
for ii=1:d
  core_kt{ii} = [ones(1, d); zeros(1, d)];
  core_kt{ii}(:, ii) = [0; 1];
  core_kt{ii}(:, mod(ii, d)+1) = [0; 1];
end
core_htd = htensor(core_kt);

disp('opts.max_rank = Inf; opts.abs_eps = 1e-5;');
disp('core_trunc = truncate(core_htd, opts);');
opts.max_rank = Inf; opts.abs_eps = 1e-5;
core_trunc = truncate(core_htd, opts);

fprintf('|| core_htd - core_trunc || = %e\n\n', ...
	norm(core_htd - core_trunc));

disp('rank(core_htd)')
rank(core_htd)
disp('rank(core_trunc)')
rank(core_trunc)

disp('core_exact = gen_exact_core(core_htd);');
core_exact = gen_exact_core(core_htd);

fprintf('|| core_htd - core_exact || = %e\n', ...
	norm(core_htd - core_exact));

disp('spy(core_exact);');
spy(core_exact);


function core_exact = gen_exact_core(core)

B = cell(1, core.nr_nodes);
U = cell(1, core.nr_nodes);
  
for ii=2:core.nr_nodes
  
  if(core.is_leaf(ii))
    U{ii} = eye(2);
  else
    ii_left  = core.children(ii, 1);
    ii_right = core.children(ii, 2);
    
    if(core.is_leaf(ii_left) && core.is_leaf(ii_right))
      
      B{ii} = dematricize(eye(4), [2 2 4], [1 2]);
      
    elseif(core.is_leaf(ii_left))
      B{ii} = zeros(2, 4, 4);
      B{ii}(1, 1, 1) = 1;
      B{ii}(2, 1, 2) = 1;
      B{ii}(1, 3, 3) = 1;
      B{ii}(2, 2, 4) = 1; B{ii}(1, 4, 4) = 1;
    elseif(core.is_leaf(ii_right))
      B{ii} = zeros(4, 2, 4);
      B{ii}(1, 1, 1) = 1;
      B{ii}(2, 1, 2) = 1;
      B{ii}(1, 2, 3) = 1;
      B{ii}(4, 1, 4) = 1; B{ii}(3, 2, 4) = 1;
    else
      B{ii} = zeros(4, 4, 4);
      B{ii}(1, 1, 1) = 1;
      B{ii}(2, 1, 2) = 1;
      B{ii}(1, 3, 3) = 1;
      B{ii}(4, 1, 4) = 1; B{ii}(3, 2, 4) = 1; B{ii}(1, 4, 4) = 1;
    end
  end
end

B{1} = zeros(4, 4);
B{1}(4, 1) = 1; B{1}(3, 2) = 1; B{1}(1, 4) = 1;
B{1}(2, 3) = 1;

core_exact = htensor(core.children, core.dim2ind, U, B);
