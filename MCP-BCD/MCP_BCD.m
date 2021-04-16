%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nonconvex penalty for tensor completion based on the MCP
%  Inputs:
%        T: a sampled tensor
%        Omega: a sampling region
%        alpha,beta,gamma,rou are parameters corresponding to MCP function
%        maxIter: the maximal iteration steps
%        epsilon: the error bound for stopping the iteration 
%  Outputs:
%        X: a reconstructed tensor
%        errlist: a list of error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, errlist] = MCP_BCD(T, Omega, alpha, beta, maxIter, epsilon, gamma, rou, X )
% The initialization of X.
if nargin < 9
    X = T;
    X(logical(1-Omega))= mean(T(Omega));
end
errlist = zeros(1,maxIter);
normT = norm(T(:));
dim = size(T);
M = cell(ndims(T), 1);
tau =alpha./beta;
betasum=sum(beta);
fprintf('MCP_BCD Iteration:     ');
for i = 1:maxIter
    for k = 1:maxIter
        fprintf('\b\b\b\b\b%5i',k);
        Xsum = 0;
        for i = 1:ndims(T)
            M{i} = Fold(MCP_minimization(Unfold(X, dim, i), tau(i),gamma,rou), dim, i);
            Xsum = Xsum + beta(i) * M{i};
        end
        Xlast = X;
        X = Xsum / betasum;
        X(Omega) = T(Omega);
        errList(k) = norm(X(:)-Xlast(:))/normT;
        if (errList(k) < 10*epsilon)
            break;
        end
    end
    error(i) = errList(k);
    if error(i) < epsilon
        errList = errList(1:k);
        break;
    end
    rou = 0.6*rou;
    beta = alpha/(0.5*rou);
    tau = alpha./beta;
    betasum = sum(beta);
    gamma = rou;
%     gamma = 0.6*gamma;
%     beta = alpha/(0.5*gamma);
%     tau = alpha./beta;
%     betasum = sum(beta);
%     rou = sqrt(2/gamma);
end
fprintf('MCP_BCD ends: total iterations = %d   difference=%f\n\n', k, errList(k));
end