function FigureBarLines(x, a)

figure;
y=0.20;
h=histogram(x,[-a:a/50:a],'Normalization','probability');
hold on;
x_sorted = sort(x);
l = size(x_sorted,1);
idx_5 = floor(0.05*l);
idx_95 = ceil(0.95*l);
x1 = x_sorted(idx_5);
x2 = x_sorted(idx_95);
% hBin=h.BinEdges(1:end-1)+h.BinWidth/2;
% Hpercent=round(h.Values/sum(h.Values)*100,2);
% W=0;
% for j1=1:100 
%    W=W+Hpercent(j1);
%    if W>=5
%       x1=j1; 
%         break
%    end
% end
% x1=a*x1/50-a;
% M=0;
% for j2=1:100 
%    M=M+Hpercent(j2);
%    if M>=95
%        x2=j2;
%         break
%    end
% end
% x2=a*(x2-50)/50;
% Hpercent2=num2cell(Hpercent);
% for i=1: length(Hpercent)
%     Hpercent2(i)={[num2str(Hpercent(i)),'%']};
% end
xlabel('($$ \bar{s} $$- s)/ s','Interpreter', 'latex')
ylabel('100%')
plot([x1,x1],[0,y],'--r','LineWidth',1.5);
plot([x2,x2],[0,y],'--r','LineWidth',1.5);
X(1)={[num2str(x1,2)]};X(2)={[num2str(x2,2)]};
text(x1*1.8,0.5*y,X(1));text(x2*1.1,0.5*y,X(2));
hold off
end