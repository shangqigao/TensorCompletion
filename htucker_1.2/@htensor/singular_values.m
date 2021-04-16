function sv = singular_values(x)
%SINGULAR_VALUES Singular values for matricizations of htensor.
%
%   SV = SINGULAR_VALUES(X) returns a cell array S containing the singular
%   values at each node in the dimension tree of an htucker X, except for
%   the root node.
%
%   See also HTENSOR, PLOT_SV, ORTHOG, GRAMIANS.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Calculate the reduced Gramians of orthogonalized x
G = gramians(orthog(x));

% Go through dimension tree
sv = cell(1, x.nr_nodes);
for ii=2:x.nr_nodes
  
  % calculate the left singular vectors U_ and singular values s
  % of X_{ii} from the gramian G{ii}.
  [tmp, sv{ii}] = htensor.left_svd_gramian(G{ii});
end