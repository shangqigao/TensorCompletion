function x_ = change_dimtree(x, children, dim2ind)
%CHANGE_DIMTREE Change dimension tree of htensor.
%
%   Y = CHANGE_DIMTREE(X, CHILDREN, DIM2IND) changes the dimension
%   tree of an htensor X to the one defined by CHILDREN, DIM2IND. If
%   this is impossible (i.e., if the dimension tree of X cannot be
%   transformed to the required one) an error message is issued. This
%   is an internal function used by CHANGE_ROOT.
%
%   See also HTENSOR, CHANGE_ROOT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Eliminate matrix at root node
ind = x.children(1, 2);

if(~x.is_leaf(ind))
  x.B{ind} = ttm(x.B{ind}, x.B{1}, 3);
else
  x.U{ind} = x.U{ind}*x.B{1}.';
end

x.B{1} = eye(size(x.B{1}, 1));

subtr = htensor.subtree(children, 1);

if(~isequal(size(children), size(x.children)) || ...
   ~isequal(sort(subtr), 1:size(x.children, 1)))
  error('CHILDREN is not consistent with X.');  
end

if(2*numel(dim2ind)-1 ~= size(children, 1) || ...
   any(any(children(dim2ind, :) ~= 0)))
  error('DIM2IND is not consistent with CHILDREN.');  
end

d = ndims(x);

% Calculate mapping between old and new indices
old2new = ones(1, 2*d-1);
new2old = ones(1, 2*d-1);

%
% x.children, x.dim2ind, old, jj
%   children,   dim2ind, new, ii
x_parent = x.parent;
for ii=size(children, 1):-1:2
  
  if(all(children(ii, :) == 0))
    jj = x.dim2ind(dim2ind == ii);    

  else  
    jj_left  = new2old(children(ii, 1));
    jj_right = new2old(children(ii, 2));
    
    par_jj_left  = x_parent(jj_left );
    par_jj_right = x_parent(jj_right);
        
    if(par_jj_left == par_jj_right)
      jj = par_jj_left;
      
    elseif(any(par_jj_left  == x.children(jj_right, :)))
      jj = par_jj_left;

    elseif(any(par_jj_right == x.children(jj_left , :)))
      jj = par_jj_right;
      
    elseif(par_jj_left == 1)
      jj = par_jj_right;

    elseif(par_jj_right == 1)
      jj = par_jj_left;
      
    else
      error('Dimension trees are not consistent.');
    end
  end
  
  old2new(jj) = ii;
  new2old(ii) = jj;
  
end

% Check mapping
if(~isequal(old2new(new2old), 1:2*d-1))
  error('Dimension trees are not consistent');
end

% Insert U, B into new cell arrays, permute as necessary
U = cell(1, 2*d-1);
B = cell(1, 2*d-1);

parent = zeros(1, size(children, 1));
ind = find(children(:, 1) ~= 0);
parent(children(ind, 1)) = ind;
parent(children(ind, 2)) = ind;

for ii=2:size(children, 1)
  jj = new2old(ii);
  
  ii1 = children(ii, 1);
  ii2 = children(ii, 2);
  
  if(ii1 == 0 && ii2 == 0)
    U{ii} = x.U{jj};
  else
    
    iipar = parent(ii);
    if(iipar == 1)
      siblings = children(iipar, :);
      iipar = siblings(siblings ~= ii);
    end
        
    jj1 = x.children(jj, 1);
    jj2 = x.children(jj, 2);
    jjpar = x.parent(jj);
    if(jjpar == 1)
      siblings = x.children(jjpar, :);
      jjpar = siblings(siblings ~= jj);
    end
    
    new_modes = new2old([ii1 ii2 iipar]);
    old_modes = [jj1 jj2 jjpar];
    
    if(old_modes(1) == new_modes(1))
      permB(1) = 1;
    elseif(old_modes(2) == new_modes(1))
      permB(1) = 2;
    else
      permB(1) = 3;
    end
    
    if(old_modes(1) == new_modes(2))
      permB(2) = 1;
    elseif(old_modes(2) == new_modes(2))
      permB(2) = 2;
    else
      permB(2) = 3;
    end
    
    if(old_modes(1) == new_modes(3))
      permB(3) = 1;
    elseif(old_modes(2) == new_modes(3))
      permB(3) = 2;
    else
      permB(3) = 3;
    end
        
    B{ii} = permute(x.B{jj}, permB);
  end
  
end

% Insert root node
ind = children(1, 1);
if(~isempty(B{ind}))
  B{1} = eye(size(B{ind}, 3));
else
  B{1} = eye(size(U{ind}, 2));
end

% Initialize new htensor
x_ = htensor(children, dim2ind, U, B, false);
