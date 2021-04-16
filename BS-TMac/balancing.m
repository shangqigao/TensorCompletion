%% The proposed balancing scheme
function [ t0,alpha ] = balancing( X,lamada )
N = ndims(X);
for i=1:N;
    x0(i)=sum(SingularValue(matricize(X,i)));
end
[x0,b0]=sort(x0);
e0=norm(x0-mean(x0))^2;
t0=cell(1,N);
for i=1:N;
    t0{i}=i;
end
t0(1:N)=t0(b0);
alpha=x0/sum(x0);
for k=1:N;
    t1=t0;
    for i=2:N;
        inset=intersect(t1{1},t1{i});
        if(isempty(inset))
            t1{1}=[t1{1} t1{i}];
            break;
        end
    end
    x1=x0;
    if x1(1)>=lamada*x1(N)
        break;
    end
    x1(1)=sum(SingularValue(matricize(X,t1{1})));
    [x1,b1]=sort(x1);
    t1(1:N)=t1(b1);
    alpha=x1/sum(x1);
%     e1=norm(x1-mean(x1))^2;
    t0=t1;
    x0=x1;
%     e0=e1;
end
end
