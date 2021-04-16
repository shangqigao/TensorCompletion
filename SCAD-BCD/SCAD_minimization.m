%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The close form solution based on the SCAD for a given matrix X
%  Inputs:
%        X: a matrix
%        lamada, gamma and rou are parameters of SCAD function
%  Output:
%        M: a solution
function [M] = SCAD_minimization(X,lamada,gamma,rou)
[U,D,V] = svd(X*X');
d = sqrt(diag(D));
g =zeros(length(d),1);
g = SCAD_thresholding(d,lamada,gamma,rou);
g(g>0) = g(g>0)./d(g>0);
M = U*diag(g)*U'*X;
end