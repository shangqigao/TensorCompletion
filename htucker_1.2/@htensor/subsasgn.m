function x = subsasgn(x, s, v)
%SUBSASGN Subscripted assignment for htensor.
%
%   X.U{I} = A assigns the matrix A to the leaf matrix at the I-th node in
%   the dimension tree of the htensor X.
%
%   X.B{I} = A assigns the 3-dimensioanl array A to the transfer tensor at
%   the I-th node.
%
%   X{I} = A assigns the array A to the leaf matrix or transfer tensor at
%   the I-th node.
%
%   The size of A must match the size of the existing leaf matrix/transfer
%   tensor. Note that the direct assignment of an element (x(2,3,1) = 3) is
%   NOT possible.
%
%   Examples:
%   x = htenrandn([5 9 7 5]); r = rank(x);
%   i = x.dim2ind(3);               % <-- get node index for 3rd dimension
%   x.U{i} = randn(7,r(i));         % <-- new leaf matrix of 3rd dimension
%   x.U{i}(6,1) = 7;                % <-- assign an individual element
%   x.B{2} = randn(r(4),r(5),r(2)); % <-- new transfer tensor
%
%   See also HTENSOR, SUBSREF, DIM2IND.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt


x.is_orthog = false;

switch s(1).type
  
 case '.'
  
  % x.U{ii} = U_new;
  if( s(1).subs == 'U' )
  
    s = s(2:end);
    str = 'U';
        
  % x.B{ii} = B_new;
  elseif( s(1).subs == 'B' )
    
    s = s(2:end);
    str = 'B';

  else
    error('Cannot change field %s, field may not exist.', s(1).subs);
  end
    
 case '()'
  error('Cannot change individual entries in an htensor.')
  
  % x{ii} = B_new; (or U_new for leaf nodes ii)
 case '{}'

   str = 'UB';
  
end

% Find index of x.U{ii}, x.B{ii} or x{ii}
ii = s(1).subs{1};

if(~isindexvector(ii) || ~isscalar(ii) || ii > x.nr_nodes)
  error('I must be one integer number between 1 and NR_NODES');
end

% Is this a leaf or non-leaf node?
is_leaf = x.is_leaf(ii);

if(strcmp(str, 'U') && ~is_leaf)
  if(~isequal(size(v), [0 0]))
    error('Index %d represents a non-leaf node', ii);
  else
    return
  end
end

if(strcmp(str, 'B') && is_leaf)
  if(~isequal(size(v), [0 0]))
    error('Index %d represents a leaf node', ii);
  else
    return
  end
end
  
s = s(2:end);

if(is_leaf)
  old_val = x.U{ii};
 else
  old_val = x.B{ii};
end

if(isempty(s))
  new_val = v;
 else
   new_val = builtin('subsasgn', old_val, s, v);
end

if( ~isequal(size(new_val), size(old_val)) )
  error('Cannot change size of leaf matrix or transfer tensor.');
end

if(is_leaf)
  x.U{ii} = new_val;
else
  x.B{ii} = new_val;
end

end

