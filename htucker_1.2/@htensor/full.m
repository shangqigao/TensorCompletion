function y = full(x)
%FULL Convert htensor to a (full) tensor.
%
%   Y = FULL(X) converts an htensor X to a (full) tensor. Use with care!
%
%   Example
%   X = htenrandn([2 4 3 5]);
%   Y = full(X) %<-- equivalent full tensor
%
%   See also HTENSOR, FULL_BLOCK.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Cell array to store the matrix U (column basis for matricization) for
% each node.
U = cell(1, x.nr_nodes);

x_is_leaf = x.is_leaf;

% Loop through all nodes, starting with the leaf nodes.
for ii=x.nr_nodes:-1:1
  
  if(x_is_leaf(ii))
    % Matrix U already known
    U{ii} = x.U{ii};
  else
    % Find child nodes
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);

    BUU = ttm(x.B{ii}, {U{ii_left}, U{ii_right}}, [1 2]);
    U{ii} = matricize(BUU, [1 2], 3, false);
    
    % Clear variables to save memory
    clear BUU;
    U{ii_left} = [];
    U{ii_right} = [];
  end 

end

% Vectorization of the full tensor is now contained in U{1}.
% Reshape to get full tensor
y = dematricize(U{1}, size(x), x.dims{1}, [], false);
