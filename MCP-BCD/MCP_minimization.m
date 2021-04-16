%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The close form solution based on the MCP for a given matrix X
%  Inputs:
%        X: a matrix
%        lamada, gamma and rou are parameters of MCP function
%  Output:
%        M: a solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M] = MCP_minimization(X,lamada,gamma,rou)
[U,D,V] = svd(X*X');
d = sqrt(diag(D));
g =zeros(length(d),1);
g = MCP_thresholding(d,lamada,gamma,rou);
g(g>0) = g(g>0)./d(g>0);
M = U*diag(g)*U'*X;
end