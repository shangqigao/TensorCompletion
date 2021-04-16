function x = full_leaves(x)
%FULL_LEAVES Convert leaf matrices U to dense matrices.
%
%  Y = FULL_LEAVES(X) results in an htensor identical to X, with
%  dense leaf matrices U{i}.
%
%  See also: SPARSE_LEAVES.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

for ii=find(x.is_leaf)
    x.U{ii} = full(x.U{ii});
end