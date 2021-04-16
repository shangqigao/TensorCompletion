function x = mrdivide(x, a)
%MRDIVIDE Scalar division for htensor.
% 
%   X = MRDIVIDE(X,A) is called when dividing an htensor X by a scalar A:
%   'X / A'.
% 
%   Example:
%   x = htenrandn([3,4,2]);
%   y = x / 5;
%
%   See also HTENSOR, HTENSOR/MTIMES

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(~isa(x, 'htensor'))
  error('First argument X must be of class htensor.');
end

if( isscalar(a) && isfloat(a) )
  x.B{1} = x.B{1}/a;
else
  error('An htensor can only be divided by a real or complex scalar.')
end
