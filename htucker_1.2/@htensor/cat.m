function x = cat(mu, x1, x2)
%CAT Concatenate two htensor objects.
%
%   X = CAT(MU, X1, X2) concatenates an htensor X1 and an htensor X2 along
%   mode MU. The tensors must have indentical sizes in all modes except in
%   MU.
%
%   Example:
%   x1 = htenrandn([3 4 7 6]); x2 = htenrandn([3 4 5 6]);
%   y = cat(3, x1, x2); % <-- returns 3x4x12x6 htensor
%
%   See also HTENSOR, CAT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin~=3 ),
  error('Requires three input arguments.');
end

% Check type of x1, x2
if( ~isa(x1, 'htensor') || ~isa(x2, 'htensor') )
  error('X1 and X2 must be of class htensor.');
end

% Check compatibility of dimension trees.
if(~equal_dimtree(x1, x2))
  error('Dimension trees of X1 and X2 differ.');
end

% Check mu.
if( ~isindexvector(mu) ),
  error('MU must be a positive integer');
end
if(mu > ndims(x1))
  error('MU exceeds tensor order.')
end

% Check sizes
sz1 = size(x1); sz1(mu) = 0;
sz2 = size(x2); sz2(mu) = 0;
if(~isequal(sz1, sz2))
  error('Tensor dimensions must agree.')
end

% Set properties of result tensor t
x = x1;
rk1 = rank(x1);
rk2 = rank(x2);
rk = rk1 + rk2;
rk(1) = 1;

x.is_orthog = false;

% Construct new leaf matrices
for ii=find(x.is_leaf)
  if(x.dim2ind(mu) ~= ii)
    x.U{ii} = [x1.U{ii}, x2.U{ii}];
  else
    x.U{ii} = blkdiag(x1.U{ii}, x2.U{ii});
  end
end

% Construct new node tensors
for ii=find(~x.is_leaf)
  
  ii_left  = x.children(ii, 1);
  ii_right = x.children(ii, 2);
  
  % Allocate space
  x.B{ii} = zeros([rk(ii_left), rk(ii_right), rk(ii)]);
  
  % Special treatment is necessary at root node
  if(ii ~= 1)
    
    x.B{ii}(1:rk1(ii_left), 1:rk1(ii_right), 1:rk1(ii)) = ...
      x1.B{ii};
    x.B{ii}(rk1(ii_left)+1:end, rk1(ii_right)+1:end, ...
	    rk1(ii)+1:end) = x2.B{ii};
  else
    
    x.B{ii}(1:rk1(ii_left)    , 1:rk1(ii_right)    ) = x1.B{ii};
    x.B{ii}(rk1(ii_left)+1:end, rk1(ii_right)+1:end) = x2.B{ii};
    
  end
end
