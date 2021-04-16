function G = gramians_cp(cp, varargin)
%GRAMIANS_CP Reduced Gramians of CP tensor.
%
%   G = GRAMIANS_CP(CP) returns a cell array containing the reduced
%   Gramians G_t for every node t of the dimension tree. CP may be a cell
%   array, or a Tensor Toolbox ktensor object. 
%   This function calls htensor to convert CP into an htensor object with
%   balanced dimension tree. When calculating the reduced Gramians, the
%   particular structure of the transfer tensors is exploited.
%
%   Note: Additional arguments can be used to specify the dimension tree,
%   as described in HTENSOR/HTENSOR.
%
%   Examples
%
%   A{1} = rand(3,3); A{2} = rand(4,3); A{3} = rand(5,3); A{4} = rand(6,3);
%   G = htensor.gramians_cp(A);      %<-- Gramians for balanced dim tree.
%   G = htensor.gramians_cp(A,'TT'); %<-- Gramians for TT-style dim tree.
%
%   See also TRUNCATE_CP, HTENSOR, GRAMIANS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(isa(cp, 'ktensor') || isa(cp, 'cell'))
  x = htensor(cp, varargin{:});
else
  error('First argument must be a cell array or ktensor object.');
end

% Child nodes
root_left  = x.children(1, 1);
root_right = x.children(1, 2);

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
    M{ii} = M{ii_left}.*M{ii_right};
    
  end
end

G = cell(1, x.nr_nodes);
G{1} = 1;

% Calculate < B{ii}, B{ii} x_1 G{ii} >_(1, 2) and _(1, 3)
G{root_left}  = M{root_right};
G{root_right}  = M{root_left};

% Traverse tree from root node downwards
for ii=find(x.is_leaf == false)
  
  if(ii == 1)
    continue;
  end
  
  % Child nodes
  ii_left  = x.children(ii, 1);
  ii_right = x.children(ii, 2);
  
  % Calculate < B{ii}, B{ii} x_1 G{ii} >_(1, 2) and _(1, 3)
  G{ii_left } = G{ii}.*M{ii_right};
  G{ii_right} = G{ii}.*M{ii_left};

end
