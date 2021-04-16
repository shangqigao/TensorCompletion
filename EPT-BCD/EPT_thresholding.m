%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The thresholding operator for the minimization based on EPT
%==========================================================================
% min_{sigmma >=0}: 1/2(sigmma-x) + lamada*f(sigmma,rou,gamma)
% where f(sigmma,rou,gamma) is the EPT function,i.e.,
% f(sigmma,rou,gamma) = 1-exp(-sigmma/gamma)
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[y] = EPT_thresholding(x,lamada,mu,gamma)
x1 = gamma*lambertw(-mu*lamada*exp(-x/gamma)/gamma^2) + x;
L_x1 = (x1-abs(x)).^2/(2*mu) + lamada*(1-exp(-x1/gamma));
L0 = x.^2/(2*mu);
y = x1;
y(x < gamma*(1+log(mu*lamada/gamma^2))) = 0;
y(L_x1 >= L0) = 0;
end
