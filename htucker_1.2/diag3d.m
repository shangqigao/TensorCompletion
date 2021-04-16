function x = diag3d(v)
%DIAG3D Return third-order diagonal tensor.
%
%   X = DIAG3D(V) returns the a third-order tensor X with diagonal elements
%   X(i, i, i) = v(i) and all other elements zero.
%   
%   Example:
%   x = diag3d([2 3]) returns
%   x(:,:,1) =
%      2     0
%      0     0
%   x(:,:,2) =
%      0     0
%      0     3
%
%   See also DIAG.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Check v
if(~isnumeric(v) || ~isvector(v))
  error('Argument v must be a (numeric) vector.');
end

% Construct n x n x n - array x
n = length(v);
x = zeros([n n n]);

% Fill in diagonal entries
dist = n^2 + n + 1;
x(1:dist:end) = v;