function [z, ind_compl_x, ind_compl_y] = ttt(x, y, ind_x, ind_y, ...
                                             ind_compl_x, ind_compl_y)
%TTT Tensor times tensor (full tensors).
%
%   TTT(X,Y) computes the outer product of two multi-dimensional arrays X
%   and Y.
% 
%   TTT(X,Y,XDIMS,YDIMS) computes the contracted product of X and Y in the
%   modes specified by the row vectors XDIMS and YDIMS. 
% 
%   TTT(X,Y,DIMS) is equivalent to calling TTT(X,Y,DIMS,DIMS).
%
%   In the case of complex tensors, the complex conjugate of X is used.
%
%   Examples
%   x = rand(4,2,3); y = rand(3,4,2);
%   z = ttt(x,y)                  %<-- outer product of x and y
%   z = ttt(x,x,1:3)              %<-- inner product of x with itself
%   z = ttt(x,y,[1 2 3],[2 3 1])  %<-- inner product of x with permuted y
%   z = ttt(x,y,[1 3],[2 1])      %<-- outer product along specified dims
% 
%   See also TTM.

%   Internal use:
%   [Z,ind_compl_x,ind_compl_y] =
%     TTT(X,Y,ind_x,ind_y,ind_compl_x,ind_compl_y)
%   with ind_x, ind_y as above (contracted modes), and ind_compl_x,
%   ind_compl_y the complement (uncontracted modes). The input arguments
%   are not checked.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt
%
% The structure of this function follows the function ttt from the Tensor
% Toolbox.

% Check input arguments
if(nargin == 6)
  % no argument checking

  % Calculate effective tensor order
  d_x = numel(ind_x) + numel(ind_compl_x);
  d_y = numel(ind_y) + numel(ind_compl_y);
  
  % Determine sizes of tensors (with trailing singleton dimensions)
  sz_x = size(x);
  sz_x(end+1:d_x) = 1;
  sz_y = size(y);
  sz_y(end+1:d_y) = 1;
else
  if(nargin < 2)
    error('Requires at least two arguments.');
  elseif(nargin == 2)
    ind_x = [];
  end
  
  if(nargin <= 3)
    ind_y = ind_x;
  end
  
  % Check x, y
  if(~isfloat(x) || ~isfloat(y))
    error('First two arguments must be multi-dimensional real/complex arrays.');
  end
  
  % Check ind_x, ind_y
  if( ~isindexvector(ind_x) || numel(unique(ind_x)) < numel(ind_x) || ...
      ~isindexvector(ind_y) || numel(unique(ind_y)) < numel(ind_y) )
    error(['IND_X and IND_Y must be vectors of positive integers,' ...
           ' without duplicate entries.']);
  end
  
  % This may help the performance of matricize; the order of the
  % dimensions in ind_x, ind_y doesn't matter otherwise
  [ind_x, sort_idx] = sort(ind_x);
  ind_y = ind_y(sort_idx);
  
  % Calculate effective tensor order
  d_x = max([ndims(x) max(ind_x)]);
  d_y = max([ndims(y), max(ind_y)]);
  
  % Determine sizes of tensors (with trailing singletons)
  sz_x = size(x);
  sz_x(end+1:d_x) = 1;
  sz_y = size(y);
  sz_y(end+1:d_y) = 1;
  
  % Compare sizes to contract
  if(~isequal(sz_x(ind_x), sz_y(ind_y)))
    error('Sizes of tensors must match for all contracted modes.');
  end
  
  % Determine dimensions that are not eliminated
  ind_compl_x = setdiff(1:d_x, ind_x);
  ind_compl_y = setdiff(1:d_y, ind_y);
  
end

% Matricize both tensors
if(any(ind_x == 1))
  transX = true;
  X = matricize(x, ind_x, ind_compl_x, false);
else
  transX = false;
  X = matricize(x, ind_compl_x, ind_x, false);
end

if(any(ind_y == 1))
  transY = false;
  Y = matricize(y, ind_y, ind_compl_y, false);
else
  transY = true;
  Y = matricize(y, ind_compl_y, ind_y, false);
end


% Calculate matricization of z, Z = X'*Y
if(transX)
  if(~transY)
    Z = X'*Y;
  else
    Z = X'*Y.';
  end
else
  if(~transY)
    Z = conj(X)*Y;
  else
    Z = conj(X)*Y.';
  end
end

% Determine size of z and dematricize
sz_z = [sz_x(ind_compl_x), sz_y(ind_compl_y)];

% Enforce at least two entries in sz_z by possibly adding singleton
% dimensions
sz_z(end+1:2) = 1;

% Dematricize matrix Z to tensor z
z = dematricize(Z, sz_z, 1:numel(ind_compl_x), ...
		numel(ind_compl_x)+1:numel(sz_z), false);
