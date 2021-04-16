function out = subsref(x, s)
%SUBSREF Subscripted reference for htensor.
%
%   Examples
%   x = htensor([5 9 3])
%   x.children returns a nr_nodes x 2 matrix specifying the dimension tree.
%   x.U returns a cell array containing all leaf bases.
%   x.B returns a cell array containing all transfer tensors.
%   x{ii} returns x.U{ii} for a leaf node ii, and x.B{ii} otherwise.
%   x(2,3,1) computes and returns the element (2,3,1) of x.
%   x(1,2:3,:) returns an htensor with the selected elements of x.
%   x(:) returns the vectorization of x.
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

switch s(1).type

 case '.'

    prop = properties('htensor');
    if( any(strcmp(s(1).subs, prop) ) ),
      out = builtin('subsref', x, s);
    else
      ll = length(prop);
      proplist = repmat({', '}, 2*ll-1, 1);
      proplist(1:2:end) = prop;
      proplist = cat(2,proplist{:});
      error(['Object htensor does not have field ' s(1).subs ...
             '. The following fields are available: ' proplist '.']);
    end
  
 case '()'

  % x(:)
  if(length(s(1).subs) == 1 && ischar(s(1).subs{1}))
    out = full(x);
    out = out(:);
  % x([1 3 4], :, 2, 1:3)
  elseif(length(s(1).subs) == ndims(x))
    
    ind = s(1).subs;
    
    for ii=find(x.is_leaf)
      % Restrict leaf matrix to the elements given by ind
      d_ii = find(x.dim2ind == ii);
      if( ~ischar(ind{d_ii}) )
        if( ~isindexvector( ind{d_ii} ) || isempty( ind{d_ii} ) )
          error('htensor only supports scalar or vector subscript indices that are positive integers.');
        end
        if( max(ind{d_ii}) > size(x.U{ii}, 1) )
          error('Subscript indices exceed htensor size.');
        end
      end
      x.U{ii} = x.U{ii}(ind{d_ii}, :);
    end    
    
    % Return scalar instead of htensor if all sizes are 1.
    if(all(size(x) == 1))
      out = full(x);
    else
      out = x;
    end
  else
    error('Number of arguments does not match order of htensor.');
  end
  
  
 case '{}'
  
  if(length(s(1).subs) ~= 1 || ~isnumeric(s(1).subs{1}) || ...
     s(1).subs{1} <= 0)
    error('{} only takes one positive integer.');
  end
  
  ii = s(1).subs{1};
  if(ii > x.nr_nodes)
    error('Index exceeds number of nodes in dimension tree.');
  end
  
  if(x.is_leaf(ii))
    out = builtin('subsref', x.U, s);
  else
    out = builtin('subsref', x.B, s);
  end
end
