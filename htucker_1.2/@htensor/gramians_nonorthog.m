function G = gramians_nonorthog(x)
%GRAMIANS_NONORTHOG Reduced Gramians of htensor.
%
%   G = GRAMIANS_NONORTHOG(X) returns a cell array containing the reduced
%   Gramians G_t for every node t of the dimension tree.
%
%   In contrast to GRAMIANS, X does not need to be (put) in orthogonalized
%   HTD.
%
%   Example
%   x = htenrandn([3 4 5 6]); G = gramians_nonorthog(x);
%
%   See also HTENSOR, TRUNCATE_NONORTHOG, GRAMIANS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Calculate M{ii} = x.U{ii}'*x.U{ii}:
M = cell(1, x.nr_nodes);

x_is_leaf = x.is_leaf;

% Traverse tree from leaf nodes upwards
for ii=x.nr_nodes:-1:2
  
  if(x_is_leaf(ii))
    % M_t = U_t' * U_t
    M{ii} = full(x.U{ii}'*x.U{ii});
  else
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);
    
    % M_t = B_t' * (M_t1 kron M_t2) * B_t
    % (interpreting B_t in matricized form)
    B_ = ttm(x.B{ii}, { M{ii_left}, M{ii_right} }, [1 2]);
    M{ii} = ttt(x.B{ii}, B_, [1 2], [1 2], 3, 3);
  end
end

G = cell(1, x.nr_nodes);
G{1} = 1;

% Traverse tree from root node downwards
for ii=find(x.is_leaf == false)
  
  % Child nodes
  ii_left  = x.children(ii, 1);
  ii_right = x.children(ii, 2);
  
  % Calculate < B{ii}, G{ii} o_1 B{ii} >_(1, 2) and _(1, 3)
  B_mod = ttm(conj(x.B{ii}), G{ii}, 3);
  
  B_mod_left  = ttm(B_mod, M{ii_right}, 2);
  B_mod_right = ttm(B_mod, M{ii_left} , 1);
  
  G{ii_left } = ttt(conj(x.B{ii}), B_mod_left, [2 3], [2 3], 1, 1);
  G{ii_right} = ttt(conj(x.B{ii}), B_mod_right, [1 3], [1 3], 2, 2);

end
