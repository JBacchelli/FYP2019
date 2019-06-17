function P=Estimate_P0(Y,X_head,PConstraints,tau_max)
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
    %%%%constraint.1, Probability constraint
       0<=P(find(PConstraints.NZIndex)>0)<=1;
    %%%%constraint.2, reachability constraint 
            P(PConstraints.NZIndex==0)==0;
    %%%%constraint.3, full observability constraint
    sum(P(:,1:length(PConstraints.OList)))==1;
    %%%%constraint.4, flow constraint 
    for jn = 1:PConstraints.OListLen
                sum( P(PConstraints.InEdge(:,jn),PConstraints.oNode_InIndex(:,jn)),1 ) - ...
                    sum( P(PConstraints.OutEdge(:,jn),PConstraints.oNode_OutIndex(:,jn)),1 ) >= 0
    end  
cvx_end
P=full(P);
P=reshape(P,[n_l,length(PConstraints.OList),tau_max]);
end

% [~,a]=size(X_head);
% [~,b]=size(Y');
% Y=Y';
% %% Estimated P by cvx
% cvx_begin  quiet
% %%optimization variable
% variable P(b,a)
% %%setting up variables relations
%     minimize(  norm(Y - X_head*P','fro' ));
%     subject to 
%     %%%%constraint.1, Probability constraint
%        0<=P(find(PConstraints.NZIndex)>0)<=1;
%     %%%%constraint.2, reachability constraint 
%             P(PConstraints.NZIndex==0)==0;
%     %%%%constraint.3, full observability constraint
%     sum(P(:,1:length(PConstraints.OList)))==1;
%     %%%%constraint.4, flow constraint 
%     for jn = 1:PConstraints.OListLen
%                 sum( P(PConstraints.InEdge(:,jn),PConstraints.oNode_InIndex(:,jn)),1 ) - ...
%                     sum( P(PConstraints.OutEdge(:,jn),PConstraints.oNode_OutIndex(:,jn)),1 ) >= 0
%     end  
% cvx_end
% P=full(P);
% P=reshape(P,[b,length(PConstraints.OList),tau_max]);
% end