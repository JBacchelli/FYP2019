function P_tilde=Toep_P(P,n_t)
[n_l,n_o,tau_max]=size(P);
HM=reshape(permute(P,[1,3,2]),[n_l*tau_max,n_o]);
P_tilde=zeros(n_l*tau_max+n_l*(n_t-1),n_o*n_t);
for tau=1:n_t
   P_tilde((tau-1)*n_l+1:(tau-1)*n_l+n_l*tau_max,(tau-1)*n_o+1:tau*n_o)=HM; 
end

% P_tilde without boundary samples
P_tilde(1:n_l*(tau_max-1),:) = 0;
P_tilde(end-n_l*(tau_max-1)+1:end,:) = 0;
end

% function P_tilde=Toep_P(P,n_t,W_P)
% [n_l,n_o,tau_max]=size(P);
% HM=reshape(permute(P,[1,3,2]),[n_l*tau_max,n_o]);
% PToep=zeros(n_l*tau_max+n_l*(n_t-1),n_o*n_t);
% 
% for tau=1:n_t
%    PToep((tau-1)*n_l+1:(tau-1)*n_l+n_l*tau_max,(tau-1)*n_o+1:tau*n_o)=HM; 
% end
% P_tilde=W_P*PToep;
% end