%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A color image for experiment: including smapling and initialization of
%parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

addpath('mylib/');
addpath('htucker_1.2');
M=double(imread('test.bmp'));
Nway = size(M);
N = ndims(M);
sr = 0.2; % sampling rate
p = round(sr*prod(Nway));
known = randsample(prod(Nway),p);
data = M(known);
Omega=zeros(Nway);
Omega(known)=1;
Omega=(Omega>0);
T=M;
T(logical(1-Omega))=255;
theta=300;
alpha=[1,1,1];
alpha = alpha / sum(alpha);
maxIter = 500;
epsilon = 1e-4;
mu = 5 * alpha ./ sqrt(size(T));
C =  0.6;
L0 = 1e-5;
coreNway=[20,20,2]; % TMac and geomCG are very sensitive to this initial rank
opts=[];
opts.tol=1e-4;
opts.alpha_adj = 0;
opts.rank_adj = -1*ones(1,3);
opts.rank_min = 5*ones(1,3);
opts.rank_max = round(1.5*coreNway);
opts.maxit=maxIter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test TMac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('TMac');
t1=tic;
[X,Y,out]=TMac(T,data,known,Nway,coreNway,opts);
time_TMac=toc(t1);
Mrec = zeros(Nway);
for i = 1:N
    Mrec = Mrec+out.alpha(i)*Fold(X{i}*Y{i},Nway,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test FaLRTC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('FaLRTC');
t2=tic;
[X_F, errList_F] = FaLRTC(...
     T,...                    % a tensor whose elements in Omega are used for estimating missing value
     Omega,...          % the index set indicating the obeserved elements
     alpha,...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
     mu,...                 % the relaxation parameters, mu >= 0. See the function for definitions
     L0,...                   % the initial step size parameter, i.e., stepsize = 1 / L;
     C,...                     % the shrinkage ratio in the range (0.5, 1). When "L" does not pass the descent enough test, we update L = L / C; 
     maxIter,...        % the maximum iterations
     epsilon,...
     data,...          % the tolerance of the relative difference of outputs of two neighbor iterations 
     known...
     );
time_FaLRTC=toc(t2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test geomCG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('geomCG');
t14=tic;
opts=[];
opts.tol=1e-4;
A = tucker_als(tensor(M), coreNway);
% create the sampling set ...
subs = makeOmegaSet( Nway, p);
% get the values of A at the sampling points ...
vals = getValsAtIndex(A, subs);
% save indices and values together in a sparse tensor
A_Omega = sptensor( subs, vals, Nway, 0);

% random initial guess:
X_init = makeRandTensor( Nway, coreNway);

[X_geomCG, ~, ~] = geomCG( A_Omega, X_init, [], opts);
X_geomCG = double(X_geomCG);
time_geomCG = toc(t14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test SpBCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('SpBCD');
t8=tic;
% beta = 1e-6*ones(1, N);
p=1/2;
beta=1e-4*ones(1,N);
[X_Sp,errList_Sp] = Sp_BCD(...
    T,...
    Omega,...
    alpha,...
    beta,...
    maxIter,...
    epsilon,...
    p...
    );
time_Sp_BCD=toc(t8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute relative error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relerr1 = norm(X_F(:)-M(:))/norm(M(:));
fprintf('FaLRTC relative error = %4.2e\n\n',relerr1);
relerr2 = norm(Mrec(:)-M(:))/norm(M(:));
fprintf('TMac relative error = %4.2e\n\n',relerr2);
relerr3 = norm(X_geomCG(:)-M(:))/norm(M(:));
fprintf('geomCG relative error = %4.2e\n\n',relerr3);
relerr4 = norm(X_Sp(:)-M(:))/norm(M(:));
fprintf('Sp-BCD relative error = %4.2e\n\n',relerr4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization the reconstructio of the color image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,3,1);
imshow(uint8(M));
title('original');
subplot(2,3,2);
imshow(uint8(T));
title('20% sampling');
subplot(2,3,3);
imshow(uint8(X_F));
title('FaLRTC');
subplot(2,3,4);
imshow(uint8(Mrec));
title('TMac');
subplot(2,3,5);
imshow(uint8(X_geomCG));
title('geomCG');
subplot(2,3,6);
imshow(uint8(X_Sp));
title('Sp-BCD');