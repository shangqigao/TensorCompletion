%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tensor Completion based on Schatten-p norm
%  Inputs:
%        T: a sampled tensor
%        Omega: a sampling region
%        alpha, beta are parameters
%        maxIter: the maximal iteration steps
%        epsilon: the error bound for stopping iteration
%        p: the paramter of Schatten-p norm
%        X: an initialization
%  Outputs:
%        X: a reconstructed tensor
%        errList: a list of error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, errList] = Sp_BCD(T, Omega, alpha,beta,maxIter, epsilon,p,X)
if nargin < 8
    X = T;
    X(logical(1-Omega)) = mean(T(Omega));
end

errList = zeros(maxIter, 1);
normT = norm(T(:));
%L = errList;
dim = size(T);
M = cell(ndims(T), 1);
tau =alpha./beta;
betasum=sum(beta);
fprintf('Sp_BCD Iteration:     ');
for k = 1:maxIter
    fprintf('\b\b\b\b\b%5i',k);
    Xsum = 0;
    for i = 1:ndims(T)
        [M{i}, ntau(i)] = Spminimization(Unfold(X, dim, i), tau(i),p);
        M{i} = Fold(M{i}, dim, i);
%         M{i} = Fold(Spminimization(Unfold(X, dim, i), tau(i),p), dim, i);
        Xsum = Xsum + beta(i) * M{i};
    end
    Xlast = X;
    X = Xsum / betasum;
    X(Omega) = T(Omega);
    errList(k) = norm(X(:)-Xlast(:))/normT;
    if (errList(k) < 10*epsilon)
        tau = ntau;
    end
    if (errList(k) < epsilon)
        errList = errList(1:k);
        break;
    end
end
fprintf('Sp_BCD ends: total iterations = %d   difference=%f\n\n', k, errList(k));
