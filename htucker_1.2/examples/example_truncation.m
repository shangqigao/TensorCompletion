function example_truncation()
%EXAMPLE_TRUNCATION Comparison of speed for different truncation methods.
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

close all;

fprintf(['This example compares the run-times of different methods of\n' ...
         'truncation. Warning! This may take a long time.\n']);
if exist('example_truncation_output.mat','file')
  fprintf('\nAlternatively, run-times from a previous run are plotted.\n\n');
  usr_str = input('Load run-times from previous run (Y/n)? ', 's');
  if(isempty(usr_str))
     usr_str = 'y';
  end
else
  usr_str = 'n';
end

figure;

if(lower(usr_str(1)) == 'y')
  load('example_truncation_output.mat');
else
  nr_summands = [2 3 4 7 9];
  
  opts.max_rank = 20;
    
  % Initialize summands
  x_cell = cell(1, 10);
  for ii=1:10
    x_cell{ii} = htenrandn(500*ones(5, 1), 'orthog', 20*ones(9, 1));
  end
  
  time_truncate_std  = zeros(1, numel(nr_summands));
  time_truncate_sum  = zeros(1, numel(nr_summands));
  time_truncate_succ = zeros(1, numel(nr_summands));
  err_truncate_std   = zeros(1, numel(nr_summands));
  err_truncate_sum   = zeros(1, numel(nr_summands));
  err_truncate_succ  = zeros(1, numel(nr_summands));
  
  for ii=numel(nr_summands):-1:1
    
    n = nr_summands(ii)
    
    x = x_cell{1};
    for jj=2:n
      x = x + x_cell{jj};
    end
    
    tic; x_ = truncate(x, opts); time_truncate_std(ii) = toc;
    tic; x_add = htensor.truncate_sum(x_cell(1:n), opts); 
    time_truncate_sum(ii) = toc;
    
    tic;
    x_succ = x_cell{1};
    for jj=2:n
      x_succ = x_succ + x_cell{jj};
      x_succ = truncate(x_succ, opts);
    end
    time_truncate_succ(ii) = toc;
    
    nrm_x = norm(orthog(x));
    err_truncate_std(ii) = norm(orthog(x - x_))/nrm_x;
    err_truncate_sum(ii) = norm(orthog(x - x_add))/nrm_x;
    err_truncate_succ(ii) = norm(orthog(x - x_succ))/nrm_x;
    
    loglog(nr_summands, time_truncate_std , 'bx-', ...
	       nr_summands, time_truncate_sum , 'rx-', ...
	       nr_summands, time_truncate_succ, 'mx-');
    
  end
end

loglog(nr_summands, time_truncate_std , 'rx-', ...
       nr_summands, time_truncate_sum , 'gx-', ...
       nr_summands, time_truncate_succ, 'bx-', ...
       'LineWidth', 2, 'MarkerSize', 8);
hold on;

loglog(nr_summands(end-1:end), 1.2*nr_summands(end-1:end).^4* ...
       time_truncate_std(end)/nr_summands(end)^4, 'r--', ...
       nr_summands(end-1:end), 1.2*nr_summands(end-1:end).^2* ...
       time_truncate_sum(end)/nr_summands(end)^2, 'g--', ...
       nr_summands(end-1:end), 1.2*nr_summands(end-1:end)* ...
       time_truncate_succ(end)/nr_summands(end), 'b--', ...
       'LineWidth', 2);
       
legend('Truncation with Alg. 6', 'Truncation with Alg. 7', ...
  'Subsequent add+truncate', ...
  'O(s^4)', 'O(s^2)', 'O(s)', 'Location', 'NorthWest');

xlabel('Number of summands');
ylabel('Runtime [sec.]');

set(gca,'YGrid','on');  set(gca,'YMinorGrid','off');

% If data is new, offer to save it into .mat file
if(lower(usr_str(1)) ~= 'y')
  usr_str = input(['Save run-times into file ' ...
		   'example_truncation_output.mat? (Y/n) '], 's');
  if(isempty(usr_str) || lower(usr_str(1)) == 'y')
    save('example_truncation_output.mat', 'nr_summands', ...
	 'time_truncate_std', 'time_truncate_sum', 'time_truncate_succ', ...
	 'err_truncate_std', 'err_truncate_sum', 'err_truncate_succ');
  end
end