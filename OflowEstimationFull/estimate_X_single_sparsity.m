function X=estimate_X_single_sparsity(Y,P0,n_t,error_bound)
%%%Estimate X by solving L1 norm minimization problem: min ||D*X||_1
%%%                                                s.t.    
%%%                                                    ||Y-P_tilde*X||_2<ErrorBound
%%%                                                       X>=0

[~,n_o]=size(P0);
%% Generate DCT tansform matrix
F=dctmtx(n_t);
%% Estimate X
cvx_begin quiet
    %% vectorize X
    variable X(n_o,n_t)
    %% main L1norm minimization
    minimize(norm(F * X',1))
    subject to
    norm(Y - P0*X, 'fro') <= error_bound;
    X >= 0;
cvx_end
end