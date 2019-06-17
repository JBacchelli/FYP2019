function [ Ao,Aod,EList,OList,ODList ] = AssignmentMatrix02( A,tau_max )
% Generate random traffic assignment matrices
% Input: 
%   A: Adjacent matrix
%   tau_max: maximum number of steps
% Output: 
%   Ao: O-Flow traffic assignment matrix
%       In Array form
%   Aod: OD flow traffic assignemnt matrix
%       In Array form

% Initialization
[~,PathM,ODList] = PathGenerate01( A,tau_max );
[PathN,~] = size(PathM);
[ODN,~] = size(ODList);
OList = union(ODList(:,1),[]);
OListLen = length(OList);

EList = AdjacentM2Edges01( A );
[EdgeN,~] = size(EList); 
Aod = zeros(EdgeN,ODN,tau_max);
Ao = zeros(EdgeN,OListLen,tau_max);

PathProb = 2 + rand(PathN,1);
PathProb = PathProb/sum(PathProb);
ODProb = zeros(ODN,1);
for odn = 1:ODN
    ODProb(odn) = sum( PathProb(ODList(odn,3):ODList(odn,4)) );
end

%% Construct OD traffic Assignment Matrix
for odn = 1:ODN
    for pn = ODList(odn,3):ODList(odn,4)
        stepN = sum(PathM(pn,:)~=0)-1;
        for sn = 1:stepN
            Edge = find( (EList(:,1)==PathM(pn,sn)) & ...
                EList(:,2)==PathM(pn,sn+1) );
            Aod(Edge,odn,sn) = Aod(Edge,odn,sn)+PathProb(pn)/ODProb(odn);
        end
    end
end

%% Construct O-Flow traffic Assignment Matrix
OProb = zeros(OListLen,1);
for oNode = 1:OListLen
    o_ODList = find( ODList(:,1)==OList(oNode) );
    OProb(oNode) = sum(ODProb(o_ODList));
    for odn = o_ODList
        for sn = 1:tau_max
            Ao(:,oNode,sn) = Ao(:,oNode,sn) + ...
                Aod(:,odn,sn)*ODProb(odn)/OProb(oNode);
        end
    end
end

end
