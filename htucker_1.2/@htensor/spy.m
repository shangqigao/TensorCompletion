function spy(x, opts)
%SPY Plot sparsity pattern of the nodes of htensor.
%
%   SPY(X) displays the dimension tree of an htensor X, with 
%   spy plots of the leaf matrices and transfer tensors.
%
%   SPY(X,OPTS) allows to choose options for the plot.
%
%   OPTS.PORTRAIT: {TRUE} | FALSE
%   Indicates orientation of the plot.
%
%   A new figure is generated, with (optional) title OPTS.TITLE.
%
%   Example:
%   x = htenrandn([4 5 6 7]);
%   spy(x+x)
%
%   See also: HTENSOR/PLOT_SV, SPY, SPY3

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end

if(nargin == 1)
  opts = struct;
elseif(~isa(opts, 'struct'))
  error('Second argument OPTS must be of class struct.');
end

% Set all fields of opts to lower case
for f=fieldnames(opts)'
opts.(lower(f{1})) = opts.(f{1});
end

if(~isfield(opts, 'portrait'))
  opts.portrait = true;
end

% New figure
if(~isfield(opts, 'title'))
  figure;
else
  figure('Name', opts.title);
end

% Calculate rectangles for each node of the tree
[pos, lvl, d_pos, d_lvl] = plot_tree(x);

% Initialize background plot:
background_axes = axes('Position', [0 0 1 1]);
set(gca,'Visible','on');

% Find non-leaf nodes
ind = find(~x.is_leaf);

% Draw lines between all nodes
if(opts.portrait)
  pos = 1 - pos;
  plot([lvl(ind), lvl(x.children(ind, 1))]', ...
       [pos(ind), pos(x.children(ind, 1))]', 'k', 'LineWidth', 2);
  hold on;
  plot([lvl(ind), lvl(x.children(ind, 2))]', ...
       [pos(ind), pos(x.children(ind, 2))]', 'k', 'LineWidth', 2);
else
  plot([pos(ind), pos(x.children(ind, 1))]', ...
       [lvl(ind), lvl(x.children(ind, 1))]', 'k', 'LineWidth', 2);
  hold on;
  plot([pos(ind), pos(x.children(ind, 2))]', ...
       [lvl(ind), lvl(x.children(ind, 2))]', 'k', 'LineWidth', 2);
end

% Set plot limits
xlim([0, 1]);
ylim([0, 1]);

% Make background gray and turn axis off
set(gca, 'Visible','off');

rotate_mode = rotate3d;

setAllowAxesRotate(rotate_mode, background_axes, false);

x_dims = x.dims;

% Loop through all nodes and spy matrices / matricizations
for ii=1:x.nr_nodes

  % Initialize axes (similar to subplot)
  if(opts.portrait)
    ax_ii = axes('OuterPosition', [lvl(ii)-d_lvl, pos(ii)-d_pos, ...
		    2*d_lvl, 2*d_pos]);
  else
    ax_ii = axes('OuterPosition', [pos(ii)-d_pos, lvl(ii)-d_lvl, ...
		    2*d_pos, 2*d_lvl]);
  end
  
  % M is the matrix U / matricization of B
  if(x.is_leaf(ii))
    % Plot matrix
    spy(x.U{ii});
    
    xlim_ = xlim; xlim([xlim_(1)+0.5, xlim_(2)-0.5]);
    ylim_ = ylim; ylim([ylim_(1)+0.5, ylim_(2)-0.5]);
    
    % Disable rotation
    setAllowAxesRotate(rotate_mode, ax_ii, false);
  else
    % Plot tensor
    spy3(x.B{ii});

    xlim_ = xlim; xlim([xlim_(1)+0.5, xlim_(2)-0.5]);
    ylim_ = ylim; ylim([ylim_(1)+0.5, ylim_(2)-0.5]);
    zlim_ = zlim; zlim([zlim_(1)+0.5, zlim_(2)-0.5]);
  end
  
  % Title: Dimensions represented by node ii
  dims_str = sprintf('%d, ', x_dims{ii});
  title(sprintf('Dim. %s', dims_str(1:end-2)), ...
	'VerticalAlignment', 'middle');
  
  % Eliminate xlabel, and ticks
  xlabel('')
  ylabel('')
  zlabel('')
  set(gca,'YTickLabel', []);
  set(gca,'ZTickLabel', []);
  set(gca,'XTickLabel', []);
  
end


% Define the visual representation of a binary tree with rectangles
% at each node: Each rectangle has center (pos(ii), lvl(ii)) and
% half sides d_pos, d_lvl. The whole tree fits into the rectangle
% [0, 1] x [0, 1].
function [pos, lvl, d_pos, d_lvl] = plot_tree(x)

% Define width 1 for each leave
width = zeros(x.nr_nodes, 1);
width(x.is_leaf) = 1;

% Add up widths for the nodes
for ii=x.nr_nodes:-1:1
    if(~x.is_leaf(ii))
        width(ii) = sum(width(x.children(ii, :)));
    end
end

% Available space for node ii
bounds = zeros(x.nr_nodes, 2);

% Position of node ii
pos = zeros(x.nr_nodes, 1);

% Root node has whole width available
bounds(1, :) = [0, width(1)];

% Go over nodes from root to leaves
for ii=1:x.nr_nodes
  
  if(~x.is_leaf(ii))
    % Child nodes
    ii_left  = x.children(ii, 1);
    ii_right = x.children(ii, 2);

    % Find position of node ii between the two child nodes' widths
    pos(ii) = bounds(ii, 1) + width(ii_left);
    
    % Calculate bounds for child nodes
    bounds(ii_left , :) = [bounds(ii, 1), pos(ii) ];
    bounds(ii_right, :) = [pos(ii), bounds(ii, 2)];
  else
    % Position is in the middle of bounds
    pos(ii) = mean(bounds(ii, :)); 
  end
end

% Put each position in the middle of its children's positions:
% Go through tree from leaves to root
for ii=x.nr_nodes:-1:1
  if(~x.is_leaf(ii))
    pos(ii) = mean(pos(x.children(ii, :)));
  end
end

% Scale to fit into [0, 1]
pos = pos/width(1);
width = width/width(1);

% Scale level to be in [0, 1]
lvl = (x.lvl' + 0.5)/(max(x.lvl)+1);

% Change order of levels
lvl = 1 - lvl;

% Half side length of each rectangle in width and level directions
d_pos = 0.5*min(width);
d_lvl = 0.5/(max(x.lvl)+1);
