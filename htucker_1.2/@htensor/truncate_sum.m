function [x, err, sv] = truncate_sum(x_cell, opts)
%TRUNCATE_SUM Truncate sum of htensor objects to lower-rank htensor.
%
%   Y = TRUNCATE(X_CELL, OPTS) truncates a sum of htensor objects in the
%   cell array X_CELL to a lower-rank htensor, according to OPTS:
%   OPTS.MAX_RANK (mandatory): maximal hierarchical rank.
%   OPTS.REL_EPS (optional): maximal relative norm-wise error.
%   OPTS.ABS_EPS (optional): maximal absolute norm-wise error.
%
%   Estimates for the overall error and for the error in each node of the
%   dimension tree are displayed if OPTS.DISP_ERRTREE is set to true.
%
%   [Y, ERR, SV] = TRUNCATE_LTR(X_CELL, OPTS) also returns the truncation
%   error and singular values in each node of the dimension tree.
%
%   Examples:
%   x_cell{1} = htenrandn([3 4 5 6]); x_cell{2} = htenrandn([3 4 5 6]);
%   opts.max_rank = 3;
%   y = truncate_sum(x_cell);
%
%   See also HTENSOR, TRUNCATE, TRUNCATE_STD, GRAMIANS_SUM, HTENSOR/PLUS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires exactly 2 arguments.')
end

% Check size(x_cell) and enforce row vector.
if(size(x_cell, 1) ~= 1 && size(x_cell, 2) == 1)
  x_cell = transpose(x_cell);
elseif(size(x_cell, 1) ~= 1)
  error('First argument must be a vector cell array.')
end

if(~iscell(x_cell) || isempty(x_cell) || ~all(cellfun( ...
      @(x)(isa(x, 'htensor')), x_cell)) )
  error('First argument must be a nonempty cell array of htensors.');
end

s = numel(x_cell);
x = x_cell{1};
sz = size(x);

% Check that size and dimension tree of all elements of x_cell are equal.
for kk = 1:s,
  if ( ~isequal(size(x_cell{kk}),sz) || ~equal_dimtree(x_cell{kk},x) )
    error('All htensor objects must have equal sizes and dimension trees.');
  end
end

if(~isa(opts, 'struct') || ~isfield(opts, 'max_rank') )
  error(['Second argument must be a MATLAB struct with at least one,' ...
	       ' field max_rank.']);
end

for ii=1:x.nr_nodes
  x.U{ii} = [];
  x.B{ii} = [];
end

% From opts.rel_eps/opts.abs_eps, calculate the error permitted for
% each truncation (there are 2d-2 truncations overall).
if(isfield(opts, 'rel_eps'))
opts.rel_eps = opts.rel_eps/sqrt(2*ndims(x) - 2);
end

if(isfield(opts, 'abs_eps'))
opts.abs_eps = opts.abs_eps/sqrt(2*ndims(x) - 2);
end

Q_cell = cell(1, x.nr_nodes);
R_cell = cell(1, x.nr_nodes);
for ii=find(x.is_leaf)
  % concatenate x_cell{:}.U{ii}
  x.U{ii} = cell2mat(cellfun(@(t)(t.U{ii}), x_cell, 'UniformOutput', false));
  
  % Calculate QR-decomposition
  [Q, R] = qr(x.U{ii}, 0);
  
  % Make sure that rank does not become zero
  if(size(R, 1) == 0)
    Q = ones(size(Q, 1), 1);
    R = ones(1, size(R, 2));
  end
  
  Q_cell{ii} = Q;
  R_cell{ii} = R;
  x.U{ii} = R;
end

% Calculate reduced Gramians of x
G = htensor.gramians_sum(x_cell);

% err represents node-wise truncation errors
err = zeros(x.nr_nodes, 1);

% Rank structure of the summands
% rank(x_cell{ii}) = blk_rank(:, ii)
blk_rank = cell2mat(cellfun(@(t)(rank(t)'), x_cell, 'UniformOutput', false));

x_is_leaf = x.is_leaf;

sv = cell(1, x.nr_nodes);

% Traverse tree from leaves to root
for ii=x.nr_nodes:-1:2

  if(x_is_leaf(ii))
    
    % Set U{ii} to Q
    x.U{ii} = Q_cell{ii};
    
  else
    % apply R{ii_left}, R{ii_right}:  
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);
    
    R_left_cell  = mat2cell(R_cell{ii_left }, ...
                            size(R_cell{ii_left }, 1), ...
                            blk_rank(ii_left , :));
    R_right_cell = mat2cell(R_cell{ii_right}, ...
                            size(R_cell{ii_right}, 1), ...
                            blk_rank(ii_right, :));
    
    x.B{ii} = cell(1, 1, s);
    for jj=1:s
      x.B{ii}{jj} = ttm(x_cell{jj}.B{ii}, {R_left_cell{jj}, ...
			                  R_right_cell{jj}}, [1 2]);
    end
    x.B{ii} = cell2mat(x.B{ii});
    
    % Matricize B{ii}
    B_mat = matricize(x.B{ii}, [1 2], 3, false);
    
    % Calculate QR decomposition
    [Q, R] = qr(B_mat, 0);
    R_cell{ii} = R;
    
    % Calculate dimensions of "tensor" Q
    tsize_new = size(x.B{ii});
    tsize_new(3) = size(Q, 2);
    
    % Reshape Q into tensor B{ii}
    x.B{ii} = dematricize(Q, tsize_new, [1 2], 3, false);
  end
  
  % Update reduced Gramian
  G{ii} = R_cell{ii}*G{ii}*R_cell{ii}';
  
  % calculate left singular vectors U_ and singular values s of X_{ii} from
  % the reduced Gramian G{ii}.
  [U_, sv{ii}] = htensor.left_svd_gramian(G{ii});
  
  % Calculate rank k to use, and expected error.
  [k, err(ii)] = htensor.trunc_rank(sv{ii}, opts);
  
  % truncate U_
  U_ = U_(:, 1:k);
  
  % Apply U_ to node ii:
  if(x_is_leaf(ii))
    x.U{ii} = x.U{ii}*U_;
  else
    x.B{ii} = ttm(x.B{ii}, U_, 3, 't');
  end
  
  R_cell{ii} = U_'*R_cell{ii};
end

% Apply R_cell to root node B
ii_left  = x.children(1, 1);
ii_right = x.children(1, 2);

R_left_cell  = mat2cell(R_cell{ii_left }, ...
                        size(R_cell{ii_left }, 1), blk_rank(ii_left , :));
R_right_cell = mat2cell(R_cell{ii_right}, ...
                        size(R_cell{ii_right}, 1), blk_rank(ii_right, :));

x.B{1} = ttm(x_cell{1}.B{1}, {R_left_cell{1}, R_right_cell{1}}, [1 2]);
for jj=2:s
  x.B{1} = x.B{1} + ttm(x_cell{jj}.B{1}, {R_left_cell{jj}, R_right_cell{jj}}, [1 2]);
end

x.is_orthog = true;

% Display the expected error in each node and overall.
if(isfield(opts, 'disp_errtree') && opts.disp_errtree == true)
  disp(x, 'truncation_error', err);
  
  % We know from theory that
  %
  % ||X - X_best|| <= ||X - X_|| <= err_bd <= factor*||X - X_best||
  %
  % and max(err) <= ||X - X_best||, therefore
  %
  % max(err_bd/factor, max(err)) <= ||X - X_best|| <= ||X - X_|| <= err_bd
  %
  % give upper and lower bounds for the best approximation as well
  % as the truncated version constructed here.
  %
  
  % Count top-level truncation only once
  err_ = err; err_(x.children(1, 1)) = 0;
  
  % Calculate upper bound and factor c from ||x - x_|| <= c ||x - x_best||
  err_bd = norm(err_); factor = sqrt(2*ndims(x)-3);
  
  fprintf(['\nLower/Upper bound for best approximation error:\n' ...
	   '%e <= ||X - X_best|| <= ||X - X_|| <= %e\n'], ...
	  max(err_bd/factor, max(err)), err_bd);
end
