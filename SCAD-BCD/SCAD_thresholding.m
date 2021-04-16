%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The thresholding operator for the minimization based on SCAD
%==========================================================================
% min_{sigmma >=0} 1/2(sigmma-x)+lamada*f(sigmma,rou,gamma)
% where f(sigmma,rou,gamma) is the MCP function,i.e.,
% f(sigmma,rou,gamma)= rou*sigmma,                           if sigmma<=rou
%                    = gamma*rou*sigmma-0.5*(sigmma^2+rou^2),if
%                                                     rou<sigmma<=gamma*rou
%                    = rou^2(gamma^2-1)/(2*(gamma-1)),  if sigmma>rou*gamma
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = SCAD_thresholding(x,lamada,gamma,rou)
L1 = (0<= x) & (x <= (lamada+1)*rou);
L2 = ((lamada+1)*rou < x) & (x <= gamma*rou);
y = x;
y(L1) = soft_thresholding(x(L1),lamada*rou);
y(L2) = soft_thresholding(x(L2),lamada*gamma*rou/(gamma-1))/(1-lamada/(gamma-1));
end
% The soft thresholding operator
function [z] = soft_thresholding(y,tau)
z = zeros(length(y),1);
z(y>tau) = y(y>tau)-tau;
z(y<-tau) = y(y<-tau)+tau; 
end