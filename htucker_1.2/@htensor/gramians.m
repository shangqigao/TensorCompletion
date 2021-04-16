function G = gramians(x)
%GRAMIANS Reduced Gramians of htensor in orthogonalized HTD.
%
%   G = GRAMIANS(X) returns a cell array containing the reduced
%   Gramians G_t for every node t of the dimension tree.
%
%   The reduced gramian satisfies X_t * X_t' = U_t * G_t * U_t', where X_t
%   is the t-matricization of X and U_t is the corresponding orthogonal
%   column basis from the HTD. Note that this method assumes that X is in
%   orthogonalized HTD. Otherwise, X is orthogonalized and a warning
%   message is issued.
%
%   Example
%   x = orthog(htenrandn([3 4 5 6])); G = gramians(x);
%
%   See also HTENSOR, GRAMIANS_NONORTHOG.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Check whether x is in orthogonalized HTD
if(~x.is_orthog)
  x = orthog(x);
  disp('htensor.gramians: Warning, htensor not in orthogonalized HTD.')
end

G = cell(1, x.nr_nodes);
G{1} = 1;

% Traverse tree from root node downwards.
for ii=find(x.is_leaf == false)

  % Child nodes
  ii_left  = x.children(ii, 1);
  ii_right = x.children(ii, 2);
  
  % Calculate contractions < B{ii}, G{ii} o_1 B{ii} >_(1, 2) and _(1, 3)
  B_mod = ttm(conj(x.B{ii}), G{ii}, 3);
  
  G{ii_left } = ttt(conj(x.B{ii}), B_mod, [2 3], [2 3], 1, 1);
  G{ii_right} = ttt(conj(x.B{ii}), B_mod, [1 3], [1 3], 2, 2);
end
