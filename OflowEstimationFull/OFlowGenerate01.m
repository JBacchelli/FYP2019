function [x]= OFlowGenerate01(tau_max_x,nn)
%%%x is of the size nOFlow-by-tau_max_x
%%%x is non-negative
%x=(1 + rand(1,nn*tau_max_x))*100;
x=rand(1,nn*tau_max_x);
%%%make non-negative
x=x-min(x);
x=reshape(x,nn,tau_max_x);
end
