function [ EList ] = AdjacentM2Edges01 (A)
%AdjacentM2Edge01: corresponding list of edges
NEdges=sum(sum(A));
EList=zeros(NEdges,2);
[EList(:,1),EList(:,2)]=find(A==1);
end