function X=estimate_X_single(Y,P0,n_t)
%%%Estimate X by solving least square problem: min ||Y-P_tilde*X||_2
%%%                                            s.t.
%%%                                                     X>=0
[~,n_o]=size(P0);
%% Estimate X by least square
    cvx_begin quiet
    variable X(n_o,n_t)
    minimize(norm(Y - P0*X, 'fro' ))
    subject to
    X>=0;
    cvx_end
end