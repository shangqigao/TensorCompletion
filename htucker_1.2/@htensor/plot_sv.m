function plot_sv(x, opts)
%PLOT_SV Plot singular value tree of htensor.
%
%   PLOT_SV(X) displays a tree containing at each node a semi-log plot of
%   the relative singular values sigma(i)/sigma(1) for the corresponding
%   matricization of the htensor X.
%
%   PLOT_SV(X,OPTS) allows to choose options for the plot.
%
%   OPTS.PORTRAIT: {TRUE} | FALSE
%   Indicates orientation of the plot.
%
%   OPTS.TICK_LABELS: TRUE | {FALSE}
%   By default, the ticks are not labeled to save space. In the x
%   direction, there are grid lines at rank 10, 20, 30, ... while they
%   appear in intervals of 10^(-2) in y direction. If there are no grid
%   lines in one direction, the first grid line is on the right / lower
%   boundary of the plot.
%
%   OPTS.STYLE: {'-b'}
%   Line style of the singular value decay plots.
%
%   A new figure is generated, with (optional) title OPTS.TITLE.
%
%   Example:
%   A = reciproc_sum(5, 100, 0.01, 1, 50); opts.max_rank = 20;
%   x = truncate(A, opts);
%   plot_sv(x);
%
%   See also: SPY, SINGULAR_VALS.

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
  fname = f{1};
  opts.(lower(fname)) = opts.(fname);
end

if(~isfield(opts, 'portrait'))
  opts.portrait = true;
end

if(~isfield(opts, 'tick_labels'))
  opts.tick_labels = false;
end

if(~isfield(opts, 'style'))
  opts.style = 'b-';
end

% Calculate singular values of x
s = singular_values(x);

% Calculate relative singular values of x
for ii=2:x.nr_nodes
    s{ii} = s{ii}/s{ii}(1);
end

% Smallest rel. singular value
min_s_ii = cell2mat(cellfun(@min, s, 'UniformOutput', false));
min_s = min([min_s_ii(min_s_ii ~= 0), 1e-2]);

% Round up max_rk to next 10^(-2*l)
log_min_s = min(2*floor(log10(min_s)/2), -2);

if(log_min_s < -20)
  log_min_s = -20;
end

% Largest rank (i.e. maximal number of singular values)
max_rk = max(cellfun('length', s));

% Round up max_rk to next multiple of 10
max_rk = 10*ceil(max_rk/10);

% New figure
if(~isfield(opts, 'title'))
  figure;
else
  figure('Name', opts.title);
end

% Calculate rectangles for each node of the tree 
[pos, lvl, d_pos, d_lvl] = plot_tree(x);

% Initialize background plot:
axes('Position', [0 0 1 1]);
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

x_dims = x.dims;

% Loop through all except root node and spy matrices / matricizations
for ii=2:x.nr_nodes
  
  % Initialize axes (similar to subplot)
  if(opts.portrait)
    axes('OuterPosition', [lvl(ii)-d_lvl, pos(ii)-d_pos, 2*d_lvl, 2*d_pos]);
  else
    axes('OuterPosition', [pos(ii)-d_pos, lvl(ii)-d_lvl, 2*d_pos, 2*d_lvl]);
  end
  
  % Plot singular values
  semilogy(s{ii}, opts.style, 'LineWidth', 2);
  
  % Title: Dimensions represented by node ii
  dims_str = sprintf('%d, ', x_dims{ii});
  title(sprintf('Dim. %s', dims_str(1:end-2)), ...
                'FontSize', 14, ...
                'VerticalAlignment', 'middle');

  axis tight;
  ylim([10^log_min_s, 1]);
  xlim([1 max_rk]);
  
  set(gca,'FontSize', 12);
  set(gca,'XTick',[1 10:10:max_rk]);
  set(gca,'YTick', 10.^(log_min_s:2:0));
  
  if(~opts.tick_labels)
    set(gca,'YTickLabel', []);
    set(gca,'XTickLabel', []);
  end
  
  set(gca,'XGrid','on')
  set(gca,'YGrid','on')
  set(gca,'XMinorGrid','off')
  set(gca,'YMinorGrid','off')
  
end

% Define the visual representation of a binary tree with rectangles
% at each node: Each rectangle has center (pos(ii), lvl(ii)) and
% half sides d_pos, d_lvl. The whole tree fits into the rectangle
% [0, 1] x [0, 1].
% The rectangle of the root node is not inside of the rectangle
function [pos, lvl, d_pos, d_lvl] = plot_tree(x)

% Define width needed at each node
width = zeros(x.nr_nodes, 1);
x_is_leaf = x.is_leaf;

% Add up widths for the nodes
for ii=x.nr_nodes:-1:1
    if(x_is_leaf(ii))
      width(ii) = 1;
    else
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
pos = pos/width(ii);
width = width/width(ii);

% Scale level to be in [0, 1]
lvl = (x.lvl' - 0.5)/max(x.lvl);

% Change order of levels
lvl = 1 - lvl;

% Adjust root node to be more clearly inside [0, 1] x [0, 1]
lvl(1) = lvl(1) - 0.6/max(x.lvl);

% Half side length of each rectangle in width and level directions
d_pos = 0.5*min(width);
d_lvl = 0.5/max(x.lvl);
