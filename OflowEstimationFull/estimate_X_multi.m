function X=estimate_X_multi(Y,P_tilde,n_o)
%%%Estimate X by solving least square problem: min ||Y-P_tilde*X||_2
%%%                                            s.t.
%%%                                                     X>=0
Y_vec = reshape(Y, numel(Y), 1);
[~,nt_no]=size(P_tilde);
n_t=nt_no/n_o;
%% Estimate X by least square
cvx_begin quiet
    variable X(nt_no)
    minimize(norm(Y_vec - P_tilde*X,2 ))
    subject to
    X>=0;
cvx_end 
    
X=reshape(X,n_o,n_t);
end