function ODFlow=OFlow2ODflow(P,X,ODList,EList)
%This function compute OD flow from input P and X
%and save the results into OD Flow.
% ODFlow: ODNum-by-tau_max_x matrix
% ODFlow(:,1):Origin nodes
% ODFlow(:,2):Destination nodes
% ODFlow(:,3):ODFlow value
ODNum=size(ODList,1);
ODFlow=zeros(ODNum,size(X,2));
P_dim=ndims(P);
if P_dim==3
    P=sum(P,3);
end
if size(P,2)~=size(X,1)
    error('Dimension of P and X do not match');
end
for odn=1:ODNum
    o_id=ODList(odn,1);
    d_id=ODList(odn,2);
    
    Edge_in=EList(:,2)==d_id;
    Edge_out=EList(:,1)==d_id;
    O_Edge_Flow=P(:,o_id)*X(o_id,:);
    Flow_in=sum(O_Edge_Flow(Edge_in,:),1);
    Flow_out=sum(O_Edge_Flow(Edge_out,:),1);
    ODFlow(odn,:)=Flow_in-Flow_out;
end
end