function FigureBar(x, left_ax, right_ax)
% myhist Codeby SimonLiang
% Email?idignew@126.com
% ??????????????

step = (right_ax - left_ax)/100;
figure;
h=histogram(x,(left_ax:step:right_ax),'Normalization','probability');
 
hold on;
 
%???????
hBin=h.BinEdges(1:end-1)+h.BinWidth/2;
%text(hBin,h.Values+max(h.Values)/25,num2cell(h.Values));
 
% %?????
% Hpercent=round(h.Values/sum(h.Values)*100,1);
%  
% %?????
% Hpercent2=num2cell(Hpercent);
% for i=1: length(Hpercent)
%     Hpercent2(i)={[num2str(Hpercent(i)),'%']};
% end
% text(hBin,h.Values+max(h.Values)/15,Hpercent2);%?????
 
%????
%title(['TotalCounts=',num2str(sum(h.Values))]);
xlabel('($$\bar{s}$$-s)/s','Interpreter', 'latex') 
hold off
end