function is_iv = isindexvector(ind)
%ISINDEXVECTOR Check whether input is index vector.
%
%   ISINDEXVECTOR(IND) is true if IND is a vector of integers larger
%   than zero.
%
%   Utility function for argument checking.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

is_iv = isnumeric(ind) && (isempty(ind) || ...
        ( isreal(ind) && isvector(ind) && all(isfinite(ind)) && ...
	  all(floor(ind) == ceil(ind)) && all(ind > 0)) );

