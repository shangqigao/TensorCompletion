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
coreNway=[50,50,2]; % TMac is very sensitive to this initial rank
opts=[];
opts.tol=1e-4;
opts.alpha_adj = 0;
opts.rank_adj = -1*ones(1,3);
opts.rank_min = 5*ones(1,3);
opts.rank_max = round(1.5*coreNway);
opts.maxit=maxIter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test S-LRTC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('S-LRTC');
t0=tic;
[X_S] = S_LRTC(...
     T,...                    % a tensor whose elements in Omega are used for estimating missing value
     Omega,...          % the index set indicating the obeserved elements   alpha,...             % the coefficient of the objective function,  i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_* 
     mu,...                 % the relaxation parameters, mu >= 0. See the function for definitions
     L0,...                   % the initial step size parameter, i.e., stepsize = 1 / L;
     C,...                     % the shrinkage ratio in the range (0.5, 1). When "L" does not pass the descent enough test, we update L = L / C; 
     maxIter,...        % the maximum iterations
     epsilon,...          % the tolerance of the relative difference of outputs of two neighbor iterations 
     coreNway,...
     alpha,...
     opts...
     );
time_S_LRTC=toc(t0);
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
%% Test NS-LRTC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t3=tic;
beta = 10*ones(1, N);
[X_NS] = NS_LRTC(...
    T,...                      % a tensor whose elements in Omega are used for estimating missing value
    Omega,...           % the index set indicating the obeserved elements
    beta,...                % the relaxation parameter. The larger, the closer to the original problem. See the function for definitions.
    maxIter,...         % the maximum iterations
    coreNway,...
    alpha,...
    opts,...
    theta...
    );
time_NS_LRTC=toc(t3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute relative error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relerr1 = norm(X_S(:)-M(:))/norm(M(:));
fprintf('S-LRTC relative error = %4.2e\n\n',relerr1);
relerr2 = norm(X_F(:)-M(:))/norm(M(:));
fprintf('FaLRTC relative error = %4.2e\n\n',relerr2);
relerr = norm(Mrec(:)-M(:))/norm(M(:));
fprintf('TMac relative error = %4.2e\n\n',relerr);
relerr3 = norm(X_NS(:)-M(:))/norm(M(:));
fprintf('NS-LRTC relative error = %4.2e\n\n',relerr3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization of the reconstructio of the color image
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
imshow(uint8(X_NS));
title('NS-LRTC');
subplot(2,3,6);
imshow(uint8(X_S));
title('S-LRTC');
