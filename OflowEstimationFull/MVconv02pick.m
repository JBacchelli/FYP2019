function Y_out=MVconv02pick(P,X)
%% Generate Y(With boundary samples)
[n_l,n_o,tau_max]=size(P);
[n_o2,n_t]=size(X);
if n_o2~=n_o
    error('dimension not match');
end
XV=reshape(X,[n_o*n_t,1]);

% Initialise complete extended P and X
P_tilde = Toep_P(P,n_t);

%% Generate Y_out (Without boundary samples)
Y_out = reshape(P_tilde * XV, n_l, n_t+tau_max-1);
end

% Original function
% function [Y_out,W_P,W_X]=MVconv02pick(P,X)
% %% Generate Y(With boundary samples)
% [m1,n1,tau_max_P]=size(P);
% [m2,tau_max_x]=size(X);
%     if m2~=n1
%     error('dimension not match');
%     end
% PM=reshape(permute(P,[1,3,2]),[m1*tau_max_P,n1]);
% XV=reshape(X,[m2*tau_max_x,1]);
% HToep=zeros(m1*tau_max_P+m1*(tau_max_x-1),m2*tau_max_x);
%     for tau=1:tau_max_x
%     HToep((tau-1)*m1+1:(tau-1)*m1+m1*tau_max_P,(tau-1)*n1+1:tau*n1)=PM;
%     end
% X_tilde=Toep_X(X,tau_max_P,1);
% P1=[];
%     for i=1:tau_max_P
%     P1=[P1,P(:,:,i)];
%     end
% P1=P1';
% Y=HToep *XV;
% Y=reshape(Y,[m1,tau_max_P+tau_max_x-1]);
% %% Generate W_X 
% Y2=X_tilde*P1;
% B_X=ones(size(Y2));
% B_X(1:tau_max_P-1,:)=0;
% B_X(tau_max_x+1:tau_max_x+tau_max_P-1,:)=0;
% W_X=diag(B_X(:,1));%%%%%generate W_x for X_tilde
% %% Generate W_P
% B_P=ones(size(Y));
% [a,b]=size(B_P);
% B_P(:,1:tau_max_P-1)=0;
% B_P(:,tau_max_x+1:tau_max_x+tau_max_P-1)=0;
% B_PVctor=[];
%     for j=1:b
%     B_PVctor=[B_PVctor;B_P(:,j)]; 
%     end
% W_P=diag(B_PVctor);%%generate W_P for P_tilde
% %% Generate Y_out (Without boundary samples)
% Y_out=B_P.*Y;%%%