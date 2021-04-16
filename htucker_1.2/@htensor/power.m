function y = power(x, p)
%POWER Element-by-element square for htensor.
%
%   Y = POWER(X, P), with p==2, returns the element-by-element square of X.
%   For all other values P, an error message is returned.
%
%   See also TIMES.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Exactly two arguments required');
end

if(~isa(x, 'htensor'))
  error('First argument must be an htensor');
end

if(p ~= 2)
  error('Only power(x, 2) (i.e x.^2) is supported');
end

y = x.*x;