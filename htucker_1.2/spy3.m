function spy3(S)
%SPY3 Plot sparsity pattern of order-3 tensor.
%
%   SPY3(S) plots the sparsity pattern of a third-order tensor S.
%
%   Example:
%   spy3(diag3d(rand(10,1)));
%
%   See also SPY, HTENSOR/SPY.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Check S
if(~isnumeric(S) || ndims(S) > 3)
  error('First argument must be a numeric array of order 3.');
end

% Calculate indices of non-zeros
[i1, i2, i3] = ind2sub(size(S), find(S));

% Plot new figure
plot3(i3, i2, i1, 'b.');

xlabel('3');
ylabel('2');
zlabel('1');

sz = size(S);
sz(end+1:3) = 1;

xlim([0, sz(3)+1]);
ylim([0, sz(2)+1]);
zlim([0, sz(1)+1]);

grid off;

set(gca, 'ydir','reverse');
set(gca, 'zdir','reverse');
set(gca, 'plotboxaspectratio',sz(end:-1:1)+1);
