function [AUniDirect, ABiDirect] = AdjacentMatrix01 (n)
% AdjacentMatrix01 for a n*n matrix, generate the adjacent matrix.
% AUnidirect : adjacent matrix for unidirectional graph
% ABidirect : adjacent matrix for bidirectional graph

At=zeros (n,n);
for nn=1:n-1
    At(nn,nn+1)=1;
end

AUniDirect=kron(At,eye(n))+kron(eye(n),At);
ABiDirect=AUniDirect+AUniDirect';
end
