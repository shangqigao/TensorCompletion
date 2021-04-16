function [ht, err, sv] = truncate(x, opts)
%TRUNCATE Truncate full tensor/htensor/CP to htensor.
%
%   Y = TRUNCATE(X, OPTS) truncates a tensor X to a lower-rank htensor
%   Y, according to OPTS:
%   OPTS.MAX_RANK (mandatory): maximal hierarchical rank.
%   OPTS.REL_EPS (optional): maximal relative norm-wise error.
%   OPTS.ABS_EPS (optional): maximal absolute norm-wise error.
%
%   This is a wrapper function that calls the appropriate truncation
%   depending on the class of X: 
%   - multidimensional array (i.e., full tensor): TRUNCATE_LTR
%   - tensor (full tensor, Tensor Toolbox):       TRUNCATE_LTR
%   - htensor (HTD decomposition):                TRUNCATE_STD
%   - ktensor (CP decomposition, Tensor Toolbox): TRUNCATE_CP
%   - cell array (CP decomposition):              TRUNCATE_CP
%   For each of these functions, further options can be specified in OPTS.
%
%   Estimates for the overall error and for the error in each node of the
%   dimension tree are displayed if OPTS.DISP_ERRTREE is set to true.
%
%   Except for an htensor input, OPTS.TREE_TYPE can always be used to
%   choose a different type of dimension tree, as described in DEFINE_TREE.
%   By default, a balanced dimension tree is chosen. Alternatively, the
%   dimension tree can be directly defined via OPTS.CHILDREN and
%   OPTS.DIM2IND.
%
%   Examples:
%   x = rand(8,8,8,8,8);
%   opts.max_rank = 3;
%   y = truncate(x,opts); %<-- htucker of rank 3
%
%   opts.max_rank = 20; opts.abs_eps = 0.01; opts.disp_errtree = true;
%   A = reciproc_sum(6, 100, 0.01, 1, 50); %<-- Function-valued CP tensor
%   x = truncate(A,opts); %<-- htucker with abs. error 0.01
%
%   See also HTENSOR, HTENSOR/TRUNCATE_LTR, HTENSOR/TRUNCATE_STD,
%   HTENSOR/TRUNCATE_CP.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires exactly 2 arguments.')
end

if(~isa(opts, 'struct') || ~isfield(opts, 'max_rank') )
  error(['Second argument must be a MATLAB struct with at least one,' ...
	       ' field max_rank.']);
end

if(isfloat(x) || isa(x, 'tensor'))
  [ht, err, sv] = htensor.truncate_ltr(x, opts);
elseif(isa(x, 'htensor'))
  [ht, err, sv] = truncate_std(x, opts);
elseif(isa(x, 'ktensor'))
  [ht, err, sv] = htensor.truncate_cp(x, opts);
elseif(isa(x, 'cell'))
  [ht, err, sv] = htensor.truncate_cp(x, opts);
else
  error('No truncation function available for %s.', class(x));
end
