%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A color image for experiment: including smapling and initialization of 
% parameters
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
coreNway=[50,50,2];
coreNwayBS=[50,50,7];
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
%% A balancing scheme for tensor completion, it's different from the original
% difiniion of the matricization of tensors. A new matricization for giving
% tensors is given.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('BS-TMac');
t0 = tic;
[t,weight] = balancing(T,0.5);
time = toc(t0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test BS-TMac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t5=tic;
[X_B,Y_B,out_B]=BS_TMac(data,known,Nway,coreNwayBS,weight,t,opts);
time_BS_TMac=toc(t5);
Mrec_BS = zeros(Nway);
for i = 1:N
    Mrec_BS = Mrec_BS+weight(i)*dematricize(X_B{i}*Y_B{i},Nway,t{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute relative error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relerr1 = norm(X_F(:)-M(:))/norm(M(:));
fprintf('FaLRTC relative error = %4.2e\n\n',relerr1);
relerr2 = norm(Mrec(:)-M(:))/norm(M(:));
fprintf('TMac relative error = %4.2e\n\n',relerr2);
relerr3 = norm(Mrec_BS(:)-M(:))/norm(M(:));
fprintf('BS-TMac relative error = %4.2e\n\n',relerr3);
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
imshow(uint8(Mrec_BS));
title('BS-TMac');
