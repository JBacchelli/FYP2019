function tdd = TDD(OD,ODTarget)
[trials,~] = size(OD);
[n_l,n_o] = size(ODTarget{1,1});

TDD = zeros(trials,1);
for i = 1:trials
    OD_i = OD{i,1};
    ODTarget_i = ODTarget{i,1};
    TDD(i,1) = abs(sum(OD_i) - sum(ODTarget_i)) / sum(ODTarget_i);
end
tdd = mean(TDD);
end