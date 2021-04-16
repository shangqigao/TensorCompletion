function[Y]=BS_SLRTC(M, Omega,mu, L, C, maxIter, epsilon,coreNway,alpha,t,opts,X1)
% initialization
if nargin < 12
    X1= M;
    X1(logical(1-Omega)) = mean(M(Omega));
end
Y=X1;
Y1 = X1;
Z1 = X1;
B = 0;

N = ndims(M);
dim = size(M);
Nway=dim;
Gx = zeros(dim);
coNway = zeros(1,N);
for n = 1:2:N
    Nway(n)=prod(dim(t{n}));
    coNway(n) = prod(dim)/Nway(n);
end
data=M(Omega);
known=find(Omega==1);

%funValueList = zeros(K+1,1);
%funValueList(1) = inf;

%%%%%%
normM = norm(M(:));
%errList(1) = norm(Y(:)-M(:)) / normM;
errList = zeros(maxIter, 1);
%%%%%%
nrmb = norm(data);
%% Parameters and defaults
if isfield(opts,'maxit')      maxit = opts.maxit;     else maxit = 500;             end
if isfield(opts,'tol')        tol = opts.tol;         else tol = 1e-4;              end
if isfield(opts,'maxT')       maxT = opts.maxT;       else maxT = 1e6;              end
%if isfield(opts,'alpha')      alpha = opts.alpha;     else alpha = ones(1,N)/N;     end
if isfield(opts,'alpha_adj') 
    alpha_adj = opts.alpha_adj;
else
    alpha_adj = 1;
end

if isfield(opts,'rank_adj') rank_adj = opts.rank_adj; else rank_adj = zeros(1,N);   end
if isfield(opts,'rank_inc') rank_inc = opts.rank_inc; else rank_inc = ones(1,N);    end
if isfield(opts,'rank_min') rank_min = opts.rank_min; else rank_min = ones(1,N);    end
if isfield(opts,'rank_max') rank_max = opts.rank_max; else rank_max = 50*ones(1,N); end
N2=sum(alpha(1:2:N));
N1=sum(alpha(2:2:N));
%% Data preprocessing and initialization
if isfield(opts,'X0')
    X2 = opts.X0;
else
    X2 = cell(1,N);
    for n = 1:2:N
        X2{n} = randn(Nway(n),coreNway(n));
    end
end

if isfield(opts,'Y0')
    Y2 = opts.Y0;
else
    Y2 = cell(1,N);
    for n = 1:2:N
        Y2{n} = randn(coreNway(n),coNway(n));
    end
end

% rescale the initial point based on number of elements
estMnrm = sqrt(nrmb^2*(prod(Nway(1:2:N))/length(known)));

% a2m = alpha.^2 ./ mu;
% ma = mu ./ alpha;
% 
% fylast = inf;

Lmax = 10*sum(1./mu);

tmp = zeros(1, N);
for i = 2:2:N
    % tmp(i) = max(SingularValue(Unfold(X, dim, i))) * alpha(i) * 0.4;
    tmp(i) = max(SingularValue(Unfold(X1, dim, i))) * alpha(i) * 0.3;%求展开后奇异值的最大值。
end
P = 1.15;
flatNum = 15;
slope = (tmp - mu) / (1-(maxIter-flatNum)^(-P));
offset = (mu*(maxIter-flatNum)^P - tmp) / ((maxIter-flatNum)^P-1); 
% offset = (mu*K^P - tmp) / ((K-flatNum)^P-1); 

mu0 = mu;
for n = 1:2:N
    X2{n} = X2{n}/norm(X2{n},'fro')*estMnrm^(Nway(n)/(Nway(n)+coNway(n)));
    Y2{n} = Y2{n}/norm(Y2{n},'fro')*estMnrm^(coNway(n)/(Nway(n)+coNway(n)));
end

X0 = X2; Y0 = Y2;

[known,id] = sort(known); data = data(id);

M2 = zeros(Nway); M2(known) = data;


sx = cell(1,N);
reschg = ones(1,N); 
reschg_tol = max(1e-2,10*tol);
rank_inc_num = sum(rank_adj==1);

res0 = zeros(1,N); TotalRes = 0;
res = res0;
for n = 1:2:N
    Mn = dematricize(X2{n}*Y2{n},dim,t{n});
    res0(n) = norm(Mn(known)-data);
    TotalRes = TotalRes+res0(n);
end
solX = ones(1,N);
Xsq = cell(1,N); Yt = cell(1,N); spI = cell(1,N);
for n = 1:2:N
    Yt{n} = Y2{n}';
end

Out.rank = coreNway;

start_time = tic;

alpha = alpha/sum(alpha);

fprintf('BS_SLRTC Iteration:     ');
for k = 1:maxIter
    fprintf('\b\b\b\b\b%5i',k);
    % update mu
    mu = max(slope / (k^P) +offset, mu0);
    % mu = mu0;
    % mu = mu0 * (K-k+1)^2
    
    a2m = alpha.^2 ./ mu;
    ma = mu ./ alpha;
    
    Ylast = Y;
    %%  test L
    %L = L*C;
    while true
        b = (1+sqrt(1+4*L*B)) / (2*L);
        X1 = b/(B+b) * Z1 + B/(B+b) * Ylast;
        % compute f'(x) namely "Gx" and f(x) namely "fx"
        Gx = Gx * 0;
        fx = 0;
        for i = 2 :2: N
            [temp, sigma2] = Truncate(matricize(X1,t{i}), ma(i));
            temp = dematricize(temp, dim, t{i});
            Gx = Gx + a2m(i) * temp;
            fx = fx + a2m(i)*(sum(sigma2) - sum(max(sqrt(sigma2)-ma(i), 0).^2));
        end
        Gx(Omega) = 0;
        % compute f(Ytest) namely fy
        Y1 = X1 - Gx / L;
        fy = 0;
        for i = 2 :2: N
            [sigma] = SingularValue(matricize(Y1,t{i}));
            fy = fy + a2m(i)*(sum(sigma.^2) - sum(max(sigma-ma(i), 0).^2));
        end
        % test if L(fx-fy) > \|Gx\|^2
        if (fx - fy)*L < sum(Gx(:).^2)
            if L > Lmax
                k
                % funValueList = funValueList(2:k);
                % disp('Exceed the Maximum Lipschitiz Constant');
                fprintf('S-LRTC: iterations = %d   difference=%f\n Exceed the Maximum Lipschitiz Constan\n\n', k, errList(k-1));
                errList = errList(1:k);
                return;
                % break;
            end
            L = L/C;
        else
             break;
        end
    end
    

    %% test Y = ? and if return
    % test whether the fylast satisfy the following condition
    % funValueList(k+1) = fy / 2;
    % Y = Ytest;
    % fylast = fy;
    % if abs(funValueList(k+1) - funValueList(k)) < epsilon
%         disp('Exceed the Minimum Tolerance');
%         break;
%     end
    
%     if fylast < fy
%         funValueList(k+1) = fylast / 2;
%     else
%         funValueList(k+1) = fy / 2;
%         Y = Ytest;
%         fylast = fy;
%         if abs(funValueList(k+1) - funValueList(k)) < epsilon
%             disp('Exceed the Minimum Tolerance');
%             break;
%         end
%     end
    
    %% update Z, Y, and B
    Z1 = Z1 - b*Gx;
    B = B+b;
%     % update (X,Y)
%     Y1=M2;
    for n = 1:2:N
        if alpha(n) > 0
            Mn = matricize(Y1,t{n});
            if solX(n)   
                X2{n} = Mn*Yt{n};
            end
            solX(n) = 1;

            Xsq{n} = X2{n}'*X2{n};
            Y2{n} = pinv(Xsq{n})*X2{n}'*Mn;
            Yt{n} = Y2{n}';

            if rank_adj(n) == -1
                rank_dec_adaptive();
            end
        end
    end
  
    % update M
%     Mn = Fold(X2{1}*Y2{1},Nway,1);
%     res(1) = norm(Mn(known)-data);
%     M2 = alpha(1)*Mn;
%     M_2=M2;
    Y=zeros(dim);
    for n = 1:2:N
        if alpha(n) > 0
        Mn = dematricize(X2{n}*Y2{n},dim,t{n});
        res(n) = norm(Mn(known)-data);
        Y = Y+(alpha(n)/N2)*Mn;
        end
    end
    % pass the true tensor M for evaluation
    if isfield(opts,'Mtr')
        Out.truerel(k) = norm(M2(:)-opts.Mtr(:))/norm(opts.Mtr(:));
    end
    
 %   M2(known) = data;
    
%     TotalRes0 = TotalRes;
%     TotalRes = 0;
%     for n = 1:2:N
%         if alpha(n) > 0
%             TotalRes = TotalRes+res(n)^2;
%         end
%     end
%     ratio = res./res0;   reschg = abs(1-ratio);    
     
    if rank_inc_num > 0
        for n = 1:2:N
            if alpha(n) > 0
            if coreNway(n) < rank_max(n) && reschg(n) < reschg_tol
                rank_inc_adaptive();
            end
            end
        end
    end
     
    % adaptively update weight
    
    if alpha_adj ~= 0
        alpha(1:2:N) =1./(res(1:2:N).^2); alpha(1:2:N) =N2*alpha(1:2:N)/sum(alpha(1:2:N));
    end
    % record how the rank estimates are updated
    Out.rank = [Out.rank; coreNway];
    
    % --- diagnostics, reporting, stopping checks ---
%     relerr1 = abs(TotalRes-TotalRes0)/(TotalRes0+1); 
%     relerr2 = sum(alpha(1:2:N).*res(1:2:N))/(N2*nrmb);
%     % reporting
%     Out.hist_rel(1,k) = relerr1;
%     Out.hist_rel(2,k) = relerr2;
    
    % check stopping criterion
%     crit = relerr1<tol;
%     if crit; nstall = nstall+1; else nstall = 0; end
%     if nstall>=3 || relerr2<tol break; end
%     if toc(start_time)>maxT; break; end;
  
%     X0 = X2; Y0 = Y2; M2 = (M2-M_2)/N2;
%     Y= zeros(Nway);
%     for i=1:2:N;
%         Y = Y+(alpha(i)/N2)*Fold(X2{i}*Y2{i},Nway,i);
%     end
%    res0 = res;
%     M2=Y;
    errList(k) =  norm(Y(known)-data) / norm(data);
    Y(known)=data;
    if k>1&&abs(errList(k)-errList(k-1)) < epsilon;
        break;
    end
%     if k>1&&abs(errList(k)-errList(k-1))<0.0001
%         break;
%     end
    %%%%%%%%%%%%
    % errList(end+1) = norm(Y(:)-M(:)) / normM;
    %%%%%%%%%%%%
end
fprintf('\n'); Out.iter = k;
Out.alpha = alpha;
% k
% funValueList = funValueList(2:k+1);

%%%%%%%%
%L = errList;
%%%%%%%%
errList = errList(1:k);
fprintf('S-LRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rank_dec_adaptive()
        % decrease the estimated rank
        sx{n} = svd(Xsq{n});
        dR = sx{n}(rank_min(n):end);
        drops = dR(1:end-1)./dR(2:end);
        [dmx,imx] = max(drops);
        rel_drp = (coreNway(n)-rank_min(n))*dmx/(sum(drops)-dmx);
        %If a large drop is found, adjust the rank estimate to imx
        if rel_drp>10
            coreNway(n) = imx+rank_min(n)-1;
            % set rank_adj(n) to 0, so only decrease the rank once
            rank_adj(n) = 0;
            [Qx,Rx] = qr(X2{n},0);
            [Qy,Ry] = qr(Y2{n}',0);
            [U,S,V] = svd(Rx*Ry');
            sigv = diag(S);
            X2{n} = Qx*U(:,1:coreNway(n))*spdiags(sigv(1:coreNway(n)),...
                0,coreNway(n),coreNway(n)); 
            X0{n} = X2{n};
            
            Yt{n} = Qy*V(:,1:coreNway(n));
            Y2{n} = Yt{n}';
            Y0{n} = Y2{n};
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rank_inc_adaptive()
        % increase the estimated rank
                   
        [Q,R] = qr(Y2{n}',0);
        for ii = 1:rank_inc(n)
            rdnx = randn(coNway(n),1);
            rdnx = rdnx-Q*(Q'*rdnx);
            rdnx = rdnx/norm(rdnx);
            Q = [Q,rdnx];
        end
        Y2{n} = Q'; Yt{n} = Q; Y0{n} = Q';
        
        coreNway(n) = coreNway(n)+rank_inc(n);
        if coreNway(n) >= rank_max(n)
            rank_inc_num = rank_inc_num - 1;
        end

        if rank_inc_num == 0
            nstall = 0;
        end

        Mn = matricize(Y1,t{n});
        X2{n} = Mn*Y2{n}';
        X0{n} = X2{n};
        solX(n) = 0;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
