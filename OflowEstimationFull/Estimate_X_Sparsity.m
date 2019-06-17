function X=Estimate_X_Sparsity(Y,P_tilde,n_o,Errorbound)
%%%Estimate X by solving L1 norm minimization problem: min ||D*X||_1
%%%                                                s.t.    
%%%                                                    ||Y-P_tilde*X||_2<ErrorBound
%%%                                                                X>=0
[~,no_nt]=size(P_tilde);
n_t=no_nt/n_o;
%% Generate DCT tansform matrix
F=dctmtx(n_t); 
N = n_o;
s = repmat('F,',1,N);
Fmtx=eval(sprintf('blkdiag(%s)',s(1:end-1)));
%% Estimate X
    cvx_begin quiet
    %% vectorilize X
    variable X(no_nt,1)
    X2=reshape(X,n_o,n_t);
    X2Vector=[];
    for i=1:n_o
    X2Vector=[X2Vector,X2(i,:)];
    end
    %% main L1norm minimization
    minimize(norm(Fmtx*X2Vector',1))
    subject to
    norm(Y(:) - P_tilde*X)<=Errorbound;
    X>=0;
    cvx_end 
    X=reshape(X,n_o,n_t);
end
