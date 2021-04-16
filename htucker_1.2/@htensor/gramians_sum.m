function G = gramians_sum(x_cell)
%GRAMIANS_SUM Reduced Gramians for sum of htensor objects.
%
%   G = GRAMIANS_SUM(X_CELL) returns a cell-array containing the reduced
%   Gramians for a sum of htensor objects in the cell array X_CELL.
%
%   Example
%   x_cell{1} = htenrandn([3 4 5 6]); x_cell{2} = htenrandn([3 4 5 6]);
%   G = htensor.gramians_sum(x_cell);
%
%   See also TRUNCATE_SUM, GRAMIANS, HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 1)
  error('Requires exactly 1 argument.')
end

if(~iscell(x_cell) || isempty(x_cell) || ~all(cellfun( ...
      @(x)(isa(x, 'htensor')), x_cell)) )
  error('X_CELL must be a nonempty cell array of htensors.');
end

% Check size(x_cell) and enforce row vector.
if(size(x_cell, 1) ~= 1 && size(x_cell, 2) == 1)
  x_cell = transpose(x_cell);
elseif(size(x_cell, 1) ~= 1)
  error('X_CELL must be a vector cell array.')
end

% Check elements of x_cell
s = numel(x_cell);
x = x_cell{1};
sz = size(x);

% Check that size and dimension tree of all elements of x_cell are equal.
for kk = 1:s,
  if ( ~isequal(size(x_cell{kk}),sz) || ~equal_dimtree(x_cell{kk},x) )
    error('All htensor objects must have equal sizes and dimension trees.');
  end
end

% Calculate M{ii} = x.U{ii}'*x.U{ii}:
M = cell(1, x.nr_nodes);

x_is_leaf = x.is_leaf;

% Traverse tree from leaf nodes upwards
for ii=x.nr_nodes:-1:2
  
  if(x_is_leaf(ii))
    % Concatenate x_cell{jj}.U{ii} for jj=1:s:
    x.U{ii} = cell2mat(cellfun(@(t)(t.U{ii}), x_cell, ...
                       'UniformOutput', false));
    
    % Block structure of x.U{ii}:
    k = cell2mat(cellfun(@(t)(size(t.U{ii}, 2)), x_cell, ...
                 'UniformOutput', false));
    
    % M{ii} is a cell array containing the blocks:
    M{ii} = mat2cell(full(x.U{ii}'*x.U{ii}), k, k);
    
  else
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);
    
    % cell array of diagonal blocks of x.B{ii}:
    B_cell = cell(s, 1);
    for jj=1:s
      B_cell{jj} = x_cell{jj}.B{ii};
    end
    
    % calculate M{ii}:
    M{ii} = cell(s, s);
    for kk=1:s
      for ll=1:s
        B1_klk = ttm(B_cell{kk}, M{ii_left}{ll, kk}, 1);
        B2_llk = ttm(B_cell{ll}, M{ii_right}{kk, ll}, 2);
        
        M{ii}{kk, ll} = ttt(B1_klk, B2_llk, [1 2], [1 2], 3, 3);
        % Reshape necessary, because 3rd dimension might be singleton,
        % which would lead to column instead of row vector
        M{ii}{kk, ll} = reshape(M{ii}{kk, ll}, ...
                                [size(B1_klk, 3), size(B2_llk, 3)]);
      end
    end
  end
end

G = cell(1, x.nr_nodes);
G{1} = 1;

% Child nodes
root_left  = x.children(1, 1);
root_right = x.children(1, 2);

% Calculate < B{ii}, B{ii} x_1 G{ii} >_(1, 2) and _(1, 3)
for kk=1:s
  for ll=1:s
    B_mod = ttm(conj(x_cell{ll}.B{1}), M{root_right}{kk, ll}, 2);
    G{root_left}{kk, ll} = ttt(conj(x_cell{kk}.B{1}), B_mod, ...
                               [2 3], [2 3], 1, 1);
    
    B_mod = ttm(conj(x_cell{ll}.B{1}), M{root_left}{kk, ll}, 1);
    G{root_right}{kk, ll} = ttt(conj(x_cell{kk}.B{1}), B_mod, ...
                                [1 3], [1 3], 2, 2);
  end
end

% Traverse tree from root node downwards
for ii=find(~x.is_leaf)
  
  if(ii == 1)
    continue;
  end
  
  % Child nodes
  ii_left  = x.children(ii, 1);
  ii_right = x.children(ii, 2);
  
  B_cell = cell(s, 1);
  for jj=1:s
    B_cell{jj} = x_cell{jj}.B{ii};
  end
  
  G{ii_left} = cell(s, s);
  for kk=1:s
    for ll=1:s
      B_mod_kkl  = ttm(conj(B_cell{kk}), G{ii}{ll, kk}, 3);
      B_right_lkl = ttm(conj(B_cell{ll}), M{ii_left }{kk, ll}, 1);
      B_left_llk  = ttm(conj(B_cell{ll}), M{ii_right}{kk, ll}, 2);
      
      G{ii_left }{kk, ll} = ttt(B_mod_kkl, B_left_llk , ...
                                [2 3], [2 3], 1, 1);
      G{ii_right}{kk, ll} = ttt(B_mod_kkl, B_right_lkl, ...
                                [1 3], [1 3], 2, 2);
    end
  end
end

for ii=2:x.nr_nodes
  G{ii} = cell2mat(G{ii});
end
