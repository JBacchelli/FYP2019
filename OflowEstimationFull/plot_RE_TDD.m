function plot_RE_TDD(OD, ODTarget, xt, xttl)
figure;
xlabel(xttl);
set(gca, 'YScale', 'linear');
set(gca, 'XScale',  'log');
points = size(OD,1);

hold on;
% xtl = cell(points,1);
% for  i=1:points
%     xtl(i) = cellstr(num2str(xt(i)));
% end
% xticks(xt);
% xticklabels(xtl);
re = zeros(points,1);
tdd = zeros(points, 1);
for k = 1:points
    re(k) = RE(OD(k), ODTarget(k));% / xt(k);
    tdd(k) = TDD(OD(k), ODTarget(k));% / xt(k);
end
yyaxis left
ylabel(strcat('RE/',xttl));
ax = get(gca, 'YAxis');
p_re = plot(xt, re, '-x', 'color', ax(1).Color, 'MarkerSize',9);
p_re.LineWidth = 1.25;
yyaxis right
ylabel(strcat('TDD/',xttl));
p_tdd = plot(xt, tdd, '-x', 'color', ax(2).Color, 'MarkerSize',9);
p_tdd.LineWidth = 1.25;

hold off;
end