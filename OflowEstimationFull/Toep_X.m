function X_tilde=Toep_X(x,tau_max)
[n_o,tau_max_x]=size(x);
XM=permute(x,[2,1]);
X_tilde=zeros(tau_max_x+tau_max-1,n_o*tau_max);

for tau=1:tau_max
   X_tilde(tau:tau-1+tau_max_x,(tau-1)*n_o+1:tau*n_o)=XM; 
end
% X_tilde without boundary samples
X_tilde(1:tau_max-1,:) = 0;
X_tilde(end-tau_max+2:end,:) = 0;
end

% function X_tilde=Toep_X(x,tau_max,W_X)
% [n_o,tau_max_x]=size(x);
% XM=permute(x,[2,1]);
% XToep=zeros(tau_max_x+tau_max-1,n_o*tau_max);
% 
% for tau=1:tau_max
%    XToep(tau:tau-1+tau_max_x,(tau-1)*n_o+1:tau*n_o)=XM; 
% end
% X_tilde=W_X*XToep;
% end