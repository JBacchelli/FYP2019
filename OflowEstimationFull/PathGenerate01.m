function [Path,PathM,ODList] = PathGenerate01(A,tau_max)
%Generate paths with maximum step number
%Note that all generated paths do not have loop
%Input:
%A:Adjacent matrix
%tau-max:maximum number of steps
%Output:
%Path: a cell
%Path {i,j}{k} gives the k-th path from node i to node j
%PathM: a matrix contains all paths, PathM(k,:) gives the kth path
%ODList: 
%ODList(k,1:2) starting and ending nodes of the k-th OD
%ODList(k,3:4) the indices of the paths in PathM for the k-th OD

%[AUniDirect, ABiDirect] = AdjacentMatrix01 (3);
%A=ABiDirect;
%tau_max=4;


[nn,~] =size(A);
Path=cell(nn,nn);


%t=1
[LList(:,1),LList(:,2)]=find(A==1);
[LLen,~]=size(LList);
for Ln=1:LLen
    Path{LList(Ln,1),LList(Ln,2)}{1}=LList(Ln,:);
end

%t>1
At1=A;
for t=2:tau_max
    At=At1*A;%
    At=At-diag(diag(At));
    clear LList;
    [LList(:,1),LList(:,2)]=find(At>0);
    [LLen,~]=size(LList);
    for Ln=1:LLen 
        oNode=LList(Ln,1);dNode=LList(Ln,2);%%%change
        Ptn=length(Path{oNode,dNode});
        K=find(At1(oNode,:)>0 & A(:,dNode)'>0);
        for k=K
          % conduct operations on all the path P_i_k^(t-1)
          Pt1N=length(Path{oNode,k});
          for Pt1n=1:Pt1N%only need the paths of length t-1
              if length(Path{oNode,k}{Pt1n})==t%check whether the path belong to this time scale
                  if sum(Path{oNode,k}{Pt1n}==dNode)>0%check weather this path passes the dnode
                      At(oNode,dNode) = At(oNode,dNode)-1;%%%%%%???
                  else
                      Ptn=Ptn+1;
                      Path{oNode,dNode}{Ptn} = [Path{oNode,k}{Pt1n},dNode];
                  end
              end
          end
        end 
    end
    At1=At;
end
%generate PathM and ODList
%initialization and memory allocation
PathNum=0;
ODNum=0;
for oNode=1:nn
    for dNode=1:nn
        odPathN=length(Path{oNode,dNode});
        if odPathN~=0
            ODNum=ODNum+1;
            PathNum=PathNum+odPathN;
        end
    end
end
PathM=zeros(PathNum,tau_max+1);
ODList=zeros(ODNum,4);

%find all paths and od pairs
Pathn=1;
ODn=1;
for oNode =1 :nn
    for dNode=1:nn
        newPathNum=length(Path{oNode,dNode});
        if newPathNum==0
            continue;
        end
        %set
        ODList(ODn,1)=oNode; ODList(ODn,2)=dNode;
        ODList(ODn,3)=Pathn;
        for newPathn=1:newPathNum
            PathM(Pathn,1:length(Path{oNode,dNode}{newPathn}))= ...
                Path{oNode,dNode}{newPathn};
            Pathn=Pathn+1;
        end
        ODList(ODn,4)=Pathn-1;
        ODn=ODn+1;
    end
end
end