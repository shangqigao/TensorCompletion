%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The close form slolution based on the EPT penalty for a given matrix X
%  Inputs:
%        X: a matrix
%        glast: an estimation of gg
%        lamada,gamma,mu are parameters
%  outputs:
%        M: a solution
%        gg is it's singular values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M,gg] = EPT_minimization(X,lamada,mu,gamma,glast)
[U,D,V] = svd(X*X');
d = sqrt(diag(D));
g =zeros(length(d),1);
g = EPT_thresholding(mu*d+(1-mu)*glast,lamada,mu,gamma);
gg = g;
g(g>0) = g(g>0)./d(g>0);
M = U*diag(g)*U'*X;
end
