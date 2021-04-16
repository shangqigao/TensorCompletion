function x = plus(x1, x2)
%PLUS Binary addition for htensor.
%
%   X = PLUS(X1, X2) adds two htensors X1 and X2 of the same size. The
%   result is an htensor of the same size but with increased hierarchical
%   ranks.
%
%   See also HTENSOR, HTENSOR/MINUS, TRUNCATE_SUM.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if( nargin~=2 )
  error('Two arguments required.');
end
% Check type of x1, x2
if( ~isa(x1, 'htensor') || ~isa(x2, 'htensor') )
  error('Both arguments must be of class htensor.');
end
% Check compatibility of dimension trees.
if(~equal_dimtree(x1, x2))
  error('Dimension trees of X1 and X2 differ.');
end
% Check sizes
if(~isequal(size(x1), size(x2)))
  error('X1 and X2 must be of identical size.')
end

% Set properties of resulting tensor x
x = x1;

r1 = rank(x1);
r2 = rank(x2);
rnew = r1 + r2;
rnew(1) = 1;

x.is_orthog = false;

% Construct new leaf matrices
for ii=find(x.is_leaf)
  x.U{ii} = [x1.U{ii}, x2.U{ii}];
end

% Construct new node tensors (except root node)
for ii=find(~x.is_leaf)
  
  ii_left  = x.children(ii, 1);
  ii_right = x.children(ii, 2);
  
  % Allocate space
  x.B{ii} = zeros([rnew(ii_left), rnew(ii_right), rnew(ii)]);
    
  if(ii ~= 1)
  x.B{ii}(1:r1(ii_left), 1:r1(ii_right), 1:r1(ii)) = x1.B{ii};
  x.B{ii}(r1(ii_left)+1:end, r1(ii_right)+1:end, r1(ii)+1:end) = x2.B{ii};
  
  else
    % Special treatment of root node
    x.B{ii}(1:r1(ii_left), 1:r1(ii_right)) = x1.B{ii};
    x.B{ii}(r1(ii_left)+1:end, r1(ii_right)+1:end) = x2.B{ii};
  end
end
