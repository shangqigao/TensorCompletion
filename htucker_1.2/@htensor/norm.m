function nrm = norm(x)
%NORM Norm of htensor.
%
%   NORM(X) returns the norm of htensor X.
%
%   If X is in orthogonalized HTD, the norm can be directly computed from
%   transfer matrix at the root node: ||B{root_node}||.
%   Otherwise, the norm is computed via the (potentially less accurate)
%   formula abs(sqrt(innerprod(X, X))).
%
%   Examples:
%   x = htenrandn([9,8,7,6]);
%   norm(x)           %<-- returns norm of x
%   norm(x-x)         %<-- should have been zero
%   norm(orthog(x-x)) %<-- a much better approximation of zero
%
%   See also HTENSOR, INNERPROD, NORM_DIFF.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% For orthogonal x, only B{1} (at the root node) is needed.
if(x.is_orthog)
  nrm = norm(x.B{1}, 'fro');
  return;
end

M = cell(x.nr_nodes, 1);

x_is_leaf = x.is_leaf;

% Traverse tree from leaves upwards.
for ii=x.nr_nodes:-1:1
  
  if(x_is_leaf(ii))
    % M_t = U_t' * U_t
    M{ii} = x.U{ii}'*x.U{ii};
  else
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);
    
    % M_t = B_t' * (M_t1 kron M_t2) * B_t
    % (interpreting B_t to be in matricized form)
    B_ = ttm(x.B{ii}, { M{ii_left}, M{ii_right} }, [1 2]);
    M{ii} = ttt(x.B{ii}, B_, [1 2], [1 2], 3, 3);
    
    % Save memory
    M{ii_left} = []; M{ii_right} = [];
  end
end

nrm = abs(sqrt(M{1}));
