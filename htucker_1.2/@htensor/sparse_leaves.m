function x = sparse_leaves(x)
%SPARSE_LEAVES Convert leaf matrices U to sparse matrices.
%
%  Y = SPARSE_LEAVES(X) results in an htensor identical to X, with sparse
%  leaf matrices U{i}.
%
%  See also: FULL_LEAVES.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

for ii=find(x.is_leaf)
  x.U{ii} = sparse(x.U{ii});
end