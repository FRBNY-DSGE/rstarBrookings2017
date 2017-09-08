function PlotStatesShaded(Time,qS)
mat_quant=squeeze(qS);
matm=mat_quant;
for i=2:size(matm,2)
    matm(:,i)=matm(:,i)-mat_quant(:,i-1);
    if(isnan(matm(end,i)))
        matm(end,i)=matm(end-1,i);
    end
end
% Generate plot
h=area(Time,matm);
r=0.8;
g=0.8;
b=0.8;
set(h,'LineStyle','none')
set(h(1),'FaceColor',[1 1 1]) % white
set(h(2),'FaceColor',.9*[r g b])
set(h(3),'FaceColor',.75*[r g b])
set(h(4),'FaceColor',.75*[r g b])
set(h(5),'FaceColor',.95*[r g b])
hold on
p1=plot(Time,qS(:,3),'--k','LineWidth',1);  % median forecast in black
hold on
plot(Time,Time*0,'k','LineWidth',.25)
hold off;
set(gcf,'Color','w')
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
axis tight; box on;
datetick('x', 17,'keeplimits', 'keepticks')
%legend([p1 p2],'Median','Realized','Location','Best')