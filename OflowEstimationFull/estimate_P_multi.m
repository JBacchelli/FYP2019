function P=estimate_P_multi(Y,X_head,PC,tau_max,n_o,sp)
%%%Estimate X by solving least square problem: min ||Y-P_tilde*X||_2
%%%                                            s.t.
%%%                                                     PConstraints 1-4
[~,no_taumax]=size(X_head);
[n_l,~]=size(Y);
Y_t=Y';
%% Estimated P by cvx
cvx_begin  quiet
%%optimization variable
variable P(n_l,no_taumax)
%%setting up variables relations
    minimize(  norm(Y_t - X_head*P','fro' ));
    subject to 
    % c4 probability 
        0<=P(PC.c4)<=1;
    % c4 reachability, speed and shortest path
        P(~PC.c4)==0;
    % c3, full observability constraint
        sum(P(:,1:n_o))==1;
    % c5 flow constraint
        for jn = 1:n_o
            sum( P(PC.c5_in_edges(:,jn),PC.c5_in_check(:,jn)),1 ) - ...
                sum( P(PC.c5_out_edges(:,jn),PC.c5_out_check(:,jn)),1 ) >= 0
        end
    % c7 duplicate counts constraint
        if sp
            for l_i = 1:n_l
                for o_i = 1:n_o
                    for tau = PC.c7_enter_link(l_i,o_i)+1:PC.c7_end_link(l_i,o_i)
                        P(l_i,o_i+tau*n_o) - ...
                            P(l_i,o_i+PC.c7_enter_link(l_i,o_i)*n_o) == 0
                    end
                end
            end
        end
    
cvx_end
P=full(P);
P=reshape(P,[n_l,n_o,tau_max]);
end