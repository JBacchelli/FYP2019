function P=estimate_P_single(Y,X0,PC,lc)
%%%Estimate X by solving least square problem: min ||Y-P_tilde*X||_2
%%%                                            s.t.
%%%
n_l = size(Y, 1);
n_o = size(X0, 1);
%n_n = size(PC.c5_inflows, 3);
%% Estimated P by cvx
cvx_begin  quiet
%%optimization variable
variable P(n_l,n_o)
%%setting up variables relations
    minimize(  norm(Y - P * X0,'fro' ));
    subject to 
    % c4 probability 
        0 <= P(PC.c4) <= lc(PC.c4);
    % c4 reachability, speed and shortest path
        P(~PC.c4) == 0;
    % c3, full observability constraint
    for o_i = 1:n_o
        sum(P(PC.c3(:,o_i),o_i) ./ lc(PC.c3(:,o_i),o_i)) == 1;
    end
    % c5 flow constraint
    %if strcmp(c5, 'old')
    for jn = 1:n_o
        sum( P(PC.c5_in_edges(:,jn),PC.c5_in_check(:,jn)) ./ lc(PC.c5_in_edges(:,jn),PC.c5_in_check(:,jn)) ) - ...
            sum( P(PC.c5_out_edges(:,jn),PC.c5_out_check(:,jn)) ./ lc(PC.c5_out_edges(:,jn),PC.c5_out_check(:,jn)) ) >= 0;
    end
    %elseif strcmp(c5, 'new')
%         for n_i = 1:n_n
%             for o_i = 1:n_o
%                 % THIS NEEDS TO BE CHANGED IF NOT SAME NUMBER OF NODES AND
%                 % ORGINS
%                 if n_i ~= o_i
%                     sum( P(PC.c5_inflows(:,o_i,n_i),o_i) ./ lc(PC.c5_inflows(:,o_i,n_i),o_i) ) - ...
%                         sum( P(PC.c5_outflows(:,o_i,n_i),o_i) ./ lc(PC.c5_outflows(:,o_i,n_i),o_i) ) >= 0;
%                 end
%             end
%         end
    %end
    
cvx_end
P=full(P);
end