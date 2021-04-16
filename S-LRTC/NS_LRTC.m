%% Non-smooth low-rank tensor completion
function[X]=NS_LRTC(M, Omega,gamma,maxIter,coreNway,alpha,opts,theta,X)
% initialization
if nargin < 9
    X= M;
    X(logical(1-Omega)) = mean(M(Omega));
end
N = ndims(M);
dim = size(M);
M1 = cell(N, 1);
Nway=dim;
coNway = zeros(1,N);
for n = 1:2:N
    coNway(n) = prod(Nway)/Nway(n);
end
data=M(Omega);
known=find(Omega==1);
%%%%%%
tau = alpha./ gamma;
errList = zeros(maxIter, 1);
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
    Mn = Fold(X2{n}*Y2{n},Nway,n);
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

fprintf('NS-LRTC Iteration:     ');
for k = 1:maxIter
    fprintf('\b\b\b\b\b%5i',k);
    for n = 1:2:N
        if alpha(n) > 0
            Mn = Unfold(X,Nway,n);
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
  
    % update M1
    Xsum = 0;
    for i = 2:2:N;
        M1{i} = Fold(Pro2TraceNorm(Unfold(X, dim, i), tau(i)), dim, i);
        Xsum = Xsum + gamma(i) * M1{i};
    end
    for i=1:2:N;
        M1{i}=Fold(X2{i}*Y2{i},Nway,i);
        Xsum=Xsum+theta*alpha(i)*M1{i};
    end
    X = Xsum / (sum(gamma(2:2:N))+theta*sum(alpha(1:2:N)));
        
    % pass the true tensor M for evaluation
    if isfield(opts,'Mtr')
        Out.truerel(k) = norm(M2(:)-opts.Mtr(:))/norm(opts.Mtr(:));
    end            
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
    errList(k) =  norm(X(known)-data) / norm(data);
    X(known)=data;
    if k>1&&abs(errList(k)-errList(k-1)) < tol;
        break;
    end

end
fprintf('\n'); Out.iter = k;
Out.alpha = alpha;
% k
% funValueList = funValueList(2:k+1);

%%%%%%%%
%L = errList;
%%%%%%%%
errList = errList(1:k);
fprintf('NS-LRTC ends: total iterations = %d   difference=%f\n\n', k, errList(k));
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

        Mn = Unfold(X,Nway,n);
        X2{n} = Mn*Y2{n}';
        X0{n} = X2{n};
        solX(n) = 0;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
