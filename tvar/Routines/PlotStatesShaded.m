function PlotStatesShadedv2(Time,qS,color,transparency)
if nargin<4
    transparency = .5;
end;
if nargin < 3
r=0.8;
g=0.8;
b=0.8;
color = [r g b];
linecolor = 'k';
else
    linecolor = color*.65;
end;

fill([Time;flip(Time)],[qS(:,1); flip(qS(:,5))],.75*color,'Linestyle','none','facealpha',transparency)
hold on
fill([Time;flip(Time)],[qS(:,2); flip(qS(:,4))],.5*color,'Linestyle','none','Linestyle','none','facealpha',transparency)
hold on
plot(Time,qS(:,3),'--','Color',linecolor,'LineWidth',1.5);
hold on
plot(Time,Time*0,'k','LineWidth',.25)
hold off;
set(gcf,'Color','w')
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
axis tight; box on;
datetick('x', 17,'keeplimits', 'keepticks')
%legend([p1 p2],'Median','Realized','Location','Best')