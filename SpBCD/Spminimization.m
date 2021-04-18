%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The closed form solution of the minimization based on Schatten-p norm
%  Inputs:
%        X: a matrix
%        tau: the regularization parameter
%        p: the parameter of Schatten-p norm
%  Outputs:
%        M: a solution
%        tau: an adaptive regularization parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M,tau]=Spminimization(X,tau,p)
[U,D,V]=svd(X*X');
d=sqrt(diag(D));
g=schatten_p(p,d,tau);
g(g>0)=g(g>0)./d(g>0);
k = length(g(g>0));
t = (p*(1-p)/2)^(1/(2-p));
if k == 0
    tau = tau*0.5;
elseif k < length(d)
    d1 = (d(k+1)/(t+p*t^(p-1)))^(2-p);
    d2 = (d(k)/(t+p*t^(p-1)))^(2-p);
    tau = 0.9*d1 + 0.1*d2;
else
    d2 = (d(k)/(t+p*t^(p-1)))^(2-p);
    tau = d2*0.1;
end
M=U*diag(g)*U'*X;
end
