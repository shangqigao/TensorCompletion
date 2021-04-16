%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nonconvex penalty for tensor completion based on the EPT
%  inputs:
%        T: a sampled tensor
%        Omega: a sampling region
%        alpha: a parameter
%        beta: a parameter
%        gamma: a parameter
%        mu: a parameter
%        maxIter: maximal iteration steps
%        epsilon: an error bound for stopping iteration
%  outputs:
%        X: a reconstructed tensor from T
%        errlist: a list of error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, errlist] = EPT_BCD(T, Omega, alpha, beta, maxIter, epsilon, mu, gamma, X )
% The initialization of X.
if nargin < 9
    X = T;
    X(logical(1-Omega))= mean(T(Omega));
end
errlist = zeros(1,maxIter);
normT = norm(T(:));
dim = size(T);
M = cell(ndims(T), 1);
glast = cell(ndims(T), 1);
for i = 1:ndims(T)
    M{i} = Unfold(X, dim, i);
    glast{i} = sqrt(svd(M{i}*M{i}'));
end
tau =alpha./beta;
betasum=sum(beta);
fprintf('EPI_BCD Iteration:     ');
for i = 1:maxIter
    for k = 1:maxIter
        fprintf('\b\b\b\b\b%5i',k);
        Xsum = 0;
        for i = 1:ndims(T)
            [M{i},glast{i}] = EPT_minimization(Unfold(X, dim, i),tau(i),mu,gamma,glast{i});
            M{i} = Fold(M{i}, dim, i);
%             M{i} = Fold(EPI_minimize(Unfold(X, dim, i),tau(i),gamma), dim, i);
            Xsum = Xsum + beta(i) * M{i};
        end
        Xlast = X;
        X = Xsum / betasum;
        X(Omega) = T(Omega);
        errList(k) = norm(X(:)-Xlast(:))/normT;
        if (errList(k) < 10*epsilon)
            errList = errList(1:k);
            break;
        end
    end
    error(i) = errList(k);
    if error(i) < epsilon
        errList = errList(1:k);
        break;
    end
    gamma = 0.4*gamma;
    beta = alpha/(10*gamma^2);
    tau = alpha./beta;
    betasum = sum(beta);
end
fprintf('EPI_BCD ends: total iterations = %d   difference=%f\n\n', k, errList(k));
end
