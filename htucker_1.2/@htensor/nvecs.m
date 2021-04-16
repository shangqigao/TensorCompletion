function U = nvecs(x, mu, r)
%NVECS Dominant left singular vectors for matricization of htensor.
%
%   U = NVECS(X, MU, R) returns a matrix containing the R leading MU-mode
%   vectors of htensor X. They correspond to the leading left singular
%   vectors of the matricization X^(MU).
%   Note that at most rank(X^(MU)) vectors are returned.
%
%   See also HTENSOR, ORTHOG, GRAMIANS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 3)
  error('Requires exactly 3 arguments.');
end

if(~isa(x, 'htensor'))
  error('First argument must be of class htensor.');
end
if(~(isindexvector(mu) && isscalar(mu)) || mu > ndims(x))
  error(['Second argument MU must be an integer not larger' ...
         ' than ndims(x).']);
end
if(~(isindexvector(r) && isscalar(r)))
  error('Third argument R must be an integer.');
end

x = orthog(x);

% Calculate the reduced Gramian of orthogonalized x
G = gramians(x);

% Leaf node for dimension mu
ii = x.dim2ind(mu);

% Calculate left singular values at node ii
U = htensor.left_svd_gramian(G{ii});
U = U(:, 1:min(r,end));
U = x.U{ii}*U;