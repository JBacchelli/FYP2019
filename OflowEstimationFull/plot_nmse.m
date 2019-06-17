function plot_nmse(lgnd, nmse_stop, varargin)
figure;
%title(ttl);
xlabel('Iteration number');
ylabel('NMSE_y');
set(gca, 'YScale', 'log');

max_iter = 0;
colors = ['r','b','g','c'];
hold on;
for k = 1:nargin-2
    nmse_history = varargin{k};
    
    if size(nmse_history,1) == 1
        nmse_history = nmse_history.nmse_history;
    end
    nmse_history = cell2mat(nmse_history);
    stop_iter = find(nmse_history<nmse_stop, 1);
    max_iter = max(max_iter, stop_iter);
    p = plot(nmse_history(1:stop_iter,1), 'color', colors(k));
    p.LineWidth = 1.5;
end
p = plot(nmse_stop*ones(max_iter,1), ':');
p.LineWidth = 2;
hold off;
legend(lgnd);
end


% y=0.20;
% h=histogram(x,[-a:a/50:a],'Normalization','probability');
% hold on;
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
% xlabel('($$ \bar{s} $$- s)/ s','Interpreter', 'latex')
% ylabel('100%')
% plot([x1,x1],[0,y],'--r','LineWidth',1.5);
% plot([x2,x2],[0,y],'--r','LineWidth',1.5);
% X(1)={[num2str(x1)]};X(2)={[num2str(x2)]};
% text(x1*1.8,0.5*y,X(1));text(x2*1.1,0.5*y,X(2));
% hold off
% end