function x = mtimes(a, b)
%MTIMES Scalar multiplication for htensor.
% 
%   C = MTIMES(A,B) is called when multiplying an htensor with a scalar.
% 
%   Example:
%   x = htenrandn([3,4,2]);
%   y = 5*x
%
%   See also HTENSOR, TTV, TTM, TTT.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Find which argument is a scalar
if( isscalar(a) )
  if( ~isfloat(a) )
    error('Scalar multiplication can only be performed with real or complex numbers.');
  end
  x = b;
  x.B{1} = a*x.B{1};
elseif( isscalar(b) )
  if( ~isfloat(b) )
    error('Scalar multiplication can only be performed with real or complex numbers.');
  end
  x = a;
  x.B{1} = b*x.B{1};
else
  error('One of the two input arguments be a scalar.')
end
