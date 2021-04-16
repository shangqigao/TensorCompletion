function disp_all(x, varargin)
%DISP_ALL Command window display of htensor.
%
%   DISP_ALL(X) displays the structure/data of an htensor object.
%
%   DISP_ALL(X,NAME) also displays the specified NAME.
%
%   DISP_ALL(X,NAME,V) also displays the entry v(ii) for every node ii in the
%   dimension tree.
%
%   See also DISP, DISPLAY, HTENSOR

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

disp(x, varargin{:});
fprintf('\n')

% ordered list of nodes for display
node_list = htensor.subtree(x.children, 1);
x_is_leaf = x.is_leaf;
x_dims = x.dims;

for ii=node_list
  dims_str = sprintf('%d ', x_dims{ii});
  dims_str = dims_str(1:end-1);
  if(x_is_leaf(ii))
    fmt = get(0,'FormatSpacing');
    format compact
    mat_str = evalc('disp(x.U{ii})');
    set(0,'FormatSpacing',fmt)
    
    fprintf('U{%d}, dims [%s]:\n%s\n', ii, dims_str, mat_str)
  else
    fmt = get(0,'FormatSpacing');
    format compact
    ten_str = evalc('disp(x.B{ii})');
    set(0,'FormatSpacing',fmt)
    
    fprintf('B{%d}, dims [%s]:\n%s\n', ii, dims_str, ten_str)
  end
end

