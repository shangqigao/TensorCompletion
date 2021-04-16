%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iteratively solving scatten-p 
function [delta] = schatten_p (p,a,lamada)
eps=1e-05;
v=(lamada*p*(1-p)/2)^(1/(2-p));
n=length(a);
delta=zeros(n,1);
for i=1:n;
    df_v=v-a(i)+lamada*p*v^(p-1);
    if df_v>=0
        delta(i)=0;
    else
        f_0=a(i)^2/2;
        x0=2*v;
        for k=1:50;
            x1=x0-(x0-a(i)+lamada*p*x0^(p-1))/(1+lamada*p*(p-1)*x0^(p-2));
            if abs(x1-x0)<eps
                break;
            end
            x0=x1;
        end
        f_x1=(x1-a(i))^2/2+lamada*x1^p;
        if f_0<f_x1
            delta(i)=0;
        else
            delta(i)=x1;
        end
    end
end
end