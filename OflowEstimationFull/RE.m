function re = RE(OD, ODTarget)

[trials,~] = size(OD);
[n_od,n_T]=size(ODTarget{1,1});
RelativeError=cell(trials,1);

max_re = 0;

for i=1:trials
    ODi=OD{i,1};
    ODTargeti=ODTarget{i,1};
    RE=zeros(n_od,n_T);
    for n=1:n_od
        for m=1:n_T
            if ODTargeti(n,m) ~= 0
                RE(n,m)=(ODi(n,m)-ODTargeti(n,m))/ODTargeti(n,m);
            end
        end
    end
    m = mean(mean(abs(RE)));
    
    RelativeError{i,1} = sqrt(1/2 * sum(sum((RE .* RE))));
end
re = mean(cell2mat(RelativeError));
end