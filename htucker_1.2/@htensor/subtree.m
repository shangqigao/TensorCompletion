function sub_ind = subtree(children, ii)
%SUBTREE Return all nodes in the subtree of a node.
%
%   T = SUBTREE(CHILDREN, I) returns a row vector containing all nodes in
%   the subtree of node I in the dimension tree. The argument CHILDREN
%   should either be an htensor object, or an N x 2 array CHILDREN defining
%   the structure of a dimension tree (as described in DEFINE_TREE).
%
%   Example
%   x = htenrandn([3 4 5 6 7 8]);
%   subtree(x,1)   %<-- returns a length-4 vector
%
%   See also: HTENSOR, DEFINE_TREE, DISP_TREE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires exactly 2 arguments.');
end

% Process argument CHILDREN.
if(isa(children, 'htensor'))
  children = children.children;
elseif( ~isnumeric(children) || size(children, 2) ~= 2 || ...
    any(floor(children(:)) ~= ceil(children(:))) || ...
    any(children(:) > size(children, 1)) || ...
    any(children(:) < 0) )
  error(['First argument CHILDREN must be an N x 2 array of integers' ...
	       ' between 0 and N, defining a binary tree.']);
end

% Check input ii
if(~(isindexvector(ii) && isscalar(ii)) || ii > size(children, 1))
  error(['Node index I must be a positive non-zero integer, and' ...
	 ' must not exceed the number of nodes.']);
end

sub_ind = ii;
ind = 1;

while(length(sub_ind) >= ind && length(sub_ind) <= size(children, 1))
  if(all(children(sub_ind(ind), :) == [0 0]))
    %do nothing
  else
    sub_ind = [sub_ind(1:ind), children(sub_ind(ind), :), sub_ind(ind+1:end)];
  end
  ind = ind + 1;
end

if(length(sub_ind) ~= length(unique(sub_ind)))
  error('The tree defined by CHILDREN contains a loop.');
end
