%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The thresholding operator for the minimization based on MCP
%==========================================================================
% min_{sigmma >=0}: 1/2(sigmma-x) + lamada*f(sigmma,rou,gamma)
% where f(sigmma,rou,gamma) is the MCP function,i.e.,
% f(sigmma,rou,gamma)= rou*sigmma-sigmma^2/(2*gamma), if sigmma<=gamma*rou
%                    = 1/2*gamma*rou^2,               if sigmma>gamma*rou 
%==========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = MCP_thresholding(x,lamada,gamma,rou)
y = x;
y(x<=rou*gamma) = soft_thresholding(x(x<=rou*gamma),lamada*rou);
end
% The soft thresholding operator
function [z] = soft_thresholding(y,tau)
z = zeros(length(y),1);
z(y>tau) = y(y>tau)-tau;
z(y<-tau) = y(y<-tau)+tau; 
end