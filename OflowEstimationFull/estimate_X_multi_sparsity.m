function X=estimate_X_multi_sparsity(Y,P_tilde,n_o,error_bound)
%%%Estimate X by solving L1 norm minimization problem: min ||D*X||_1
%%%                                                s.t.    
%%%                                                    ||Y-P_tilde*X||_2<ErrorBound
%%%                                                                X>=0
Y_vec = reshape(Y, numel(Y), 1);
[~,nt_no]=size(P_tilde);
n_t=nt_no/n_o;
%% Generate DCT tansform matrix
F=dctmtx(n_t);
s = repmat('F,',1,n_o);
Fmtx=eval(sprintf('blkdiag(%s)',s(1:end-1)));
%% Estimate X
cvx_begin quiet
    %% vectorilize X
    variable X(nt_no,1)
    X2=reshape(X,n_o,n_t);
    X2Vector = cvx(zeros(1,nt_no));
    for i=1:n_o
        X2Vector(1+n_t*(i-1):i*n_t)=X2(i,:);
    end
    %% main L1norm minimization
    minimize(norm(Fmtx*X2Vector',1))
    subject to
    norm(Y(:) - P_tilde*X) <= error_bound;
    X >= 0;
cvx_end 
X = reshape(X,n_o,n_t);
end
% %% Estimate X
% cvx_begin quiet
%     %% vectorilize X
%     variable X(nt_no)
%     %% main L1 norm minimization
%     minimize(norm(Fmtx*X,1))
%     subject to
%     norm(Y_vec(:) - P_tilde*X) <= error_bound;
%     X >= 0;
% cvx_end
% X = reshape(X,n_o,n_t);
% end