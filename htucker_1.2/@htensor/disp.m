function disp(x, name, v)
%DISP Command window display of dimension tree of htensor.
%
%   DISP(X) displays the structure of an htensor object.
%
%   DISP(X,NAME) also displays the specified NAME.
%
%   DISP(X,NAME,V) also displays the entry v(ii) for every node ii in
%   the dimension tree.
%
%   See also DISP_ALL, DISPLAY, HTENSOR

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

if (nargin < 2 || ~ischar(name))
  name = 'ans';
end

if(nargin == 3 && (~isvector(v) || numel(v) ~= x.nr_nodes))
  error('V must be a vector with X.NR_NODES elements.');
end

sz = size(x);
size_str = [sprintf('%d x ', sz(1:end-1)), sprintf('%d', sz(end))];
fprintf(1, '%s is an htensor of size %s\n', name, size_str);

% calculate maximum number of digits
d = ndims(x);
dims_digits = ceil(log10(d));
node_digits = ceil(log10(2*d-1));
rank_digits = ceil(log10(max(rank(x))));

f_node = ['%' num2str(node_digits) 'd'];
f_rank = ['%' num2str(rank_digits) 'd'];

% ordered list of nodes for display
node_list = htensor.subtree(x.children, 1);

% Precompute x.lvl
lvl = x.lvl;

% loop over nodes
for ii=node_list
  
  % calculate indent according to level
  dims = x.dims{ii};
  indent = repmat('  ', 1, lvl(ii));
  
  % string indicating the dimensions of the node
  if(length(dims) > 1)
    out_dims = sprintf('%s%d-%d', indent, min(dims), max(dims));
  else
    out_dims = sprintf('%s%d', indent, dims);
  end
  
  % balance the indent, to get a second column
  max_indent = 2*max(lvl) + 2*dims_digits+1;
  indent_balance = repmat(' ', 1, 3 + max_indent - length(out_dims));
  
  % string with node number and size of U / B
  if(length(dims) == 1)
    out_size = sprintf([f_node '; ' f_rank ' ' f_rank], ...
                       ii, size(x.U{ii}, 1), size(x.U{ii}, 2));
  else
    out_size = sprintf([f_node '; ' f_rank ' ' f_rank ' ' f_rank], ...
                       ii, size(x.B{ii}, 1), ...
                       size(x.B{ii}, 2), size(x.B{ii}, 3));
  end

  % balance the indent, to get a third column
  max_indent = (node_digits+2) + 3*(rank_digits+2) + 1;
  indent_balance2 = repmat(' ', 1, 3 + max_indent - length(out_size));
  
  if(nargin <= 2)
    out_v = '';
  else
    out_v = sprintf('%d', v(ii));
  end
  
  fprintf(1, '%s%s%s%s%s\n', out_dims, indent_balance, out_size, ...
	        indent_balance2, out_v);
end
