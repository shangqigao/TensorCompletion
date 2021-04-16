function x = uminus(x)
%UMINUS Unary minus (-) of htensor.
%
%   See also HTENSOR.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Change sign of tranfer matrix at root node
x.B{1} = -x.B{1};