function PConstraints = GetPConstraintsMultiStep21( A,tau_max,P_t )
% This is for multi step rigid model
%   A is the Adjacent Matrix
%   tau_max is the number of steps inherently in the assignment matrix P_t
%   P_t is an array for assignment matrices, where the 3rd dimension is
%      time. 

% Initialization
[~,~,ODList] = PathGenerate01( A,tau_max );
OList = union(ODList(:,1),[]);
OListLen = length(OList);

EList = AdjacentM2Edges01( A );
[EdgeN,~] = size(EList); 

if ndims(P_t)<3
    error('Dimension of P_t must be 3');
elseif size(P_t,3)~=tau_max
    error('The 3rd dimension of P_t must be tau_max');
elseif tau_max <= 1
    error('tau_max should be at least 2 for multi-step models');
end
[P_m,P_n,~] = size(P_t);
P_matrix = reshape( P_t, [P_m,P_n*tau_max] );
% All the constraints are defined w.r.t. P_matrix

% Probability constraint & Speed constraint
PConstraints.NZIndex = P_matrix~=0;

% Observability Constraints
PConstraints.OList = OList;
PConstraints.OListLen = OListLen;
PConstraints.StartingEdge = repmat(1==0,P_m,P_n);
for oNn = 1:OListLen
    oNode = OList(oNn);
    PConstraints.StartingEdge(1:EdgeN,oNn) = EList(1:EdgeN,1)==oNode; 
end

%% Flow Constraints
% for each oNode, we need to find a list of middle nodes
% each middle node will corresponding to one inequality constraint where
% the first summ is for in-flow and the 2nd summ is for out-flow

PConstraints.j_InIndex = repmat(1==0,OListLen*tau_max,OListLen);
PConstraints.j_OutIndex = repmat(1==0,OListLen*tau_max,OListLen);
PConstraints.InEdge = repmat(1==0,EdgeN,OListLen);
PConstraints.OutEdge = repmat(1==0,EdgeN,OListLen);
for j_p_n = 1:OListLen
    j_p = OList(j_p_n);
    
    PConstraints.InEdge(:,j_p_n) = EList(:,2)==j_p;
    PConstraints.OutEdge(:,j_p_n) = EList(:,1)==j_p;
    
    PConstraints.oNode_InIndex(:,j_p_n) = ... 
        [ repmat(OList~=j_p,tau_max-1,1) ; repmat(1==0,OListLen,1) ];
    PConstraints.oNode_OutIndex(:,j_p_n) = ... 
        [ repmat(1==0,OListLen,1) ; repmat(OList~=j_p,tau_max-1,1) ]; 
end

end
