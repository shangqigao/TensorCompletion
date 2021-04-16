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
coreNway=[20,20,2]; % TMac is very sensitive to this initial rank
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
%% Test EPT-BCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('EPT-BCD');
t10=tic;
mu=0.9;
gamma=1e3; %for color image
% gamma=1e4; %for MRI
beta=alpha/(10*gamma^2);
[X_EPT, errlist_EPT] = EPT_BCD(...
    T,...
    Omega,...
    alpha,...
    beta,...
    maxIter,...
    epsilon,...
    mu,...
    gamma...
    );
time_EPT=toc(t10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test MCP-BCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('MCP-BCD');
t11=tic;
rou=100; %for color image
% rou=200;  %for cardiac
% rou=300;  % for video, knee
beta = alpha/(0.5*rou);  
gamma=rou;
[X_MCP, errlist_MCP] = MCP_BCD(...
    T,...
    Omega,...
    alpha,...
    beta,...
    maxIter,...
    epsilon,...
    gamma,...
    rou...
    );
time_MCP=toc(t11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test SCAD-BCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('SCAD-BCD');
t12=tic;
rou=100; %for color image
% rou=1000;  %for brain
% rou=200;  %for cardiac
% rou=300;  %for video, knee
beta = alpha/(0.5*rou);
gamma=rou;
[X_SCAD, errlist_SCAD] = SCAD_BCD(...
    T,...
    Omega,...
    alpha,...
    beta,...
    maxIter,...
    epsilon,...
    gamma,...
    rou...
    );
time_SCAD=toc(t12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute relative error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relerr1 = norm(X_F(:)-M(:))/norm(M(:));
fprintf('FaLRTC relative error = %4.2e\n\n',relerr1);
relerr2 = norm(Mrec(:)-M(:))/norm(M(:));
fprintf('TMac relative error = %4.2e\n\n',relerr2);
relerr3 = norm(X_MCP(:)-M(:))/norm(M(:));
fprintf('MCP_BCD relative error = %4.2e\n\n',relerr3);
relerr4 = norm(X_SCAD(:)-M(:))/norm(M(:));
fprintf('SCAD_BCD relative error = %4.2e\n\n',relerr4);
relerr5 = norm(X_EPT(:)-M(:))/norm(M(:));
fprintf('EPT_BCD relative error = %4.2e\n\n',relerr5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization of the reconstruction of the color image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,4,1);
imshow(uint8(M));
title('original');
subplot(2,4,2);
imshow(uint8(T));
title('20% sampling');
subplot(2,4,3);
imshow(uint8(X_F));
title('FaLRTC');
subplot(2,4,4);
imshow(uint8(Mrec));
title('TMac');
subplot(2,4,5);
imshow(uint8(X_EPT));
title('EPT-BCD');
subplot(2,4,6);
imshow(uint8(X_MCP));
title('MCP-BCD')
subplot(2,4,7);
imshow(uint8(X_SCAD));
title('SCAD-BCD');