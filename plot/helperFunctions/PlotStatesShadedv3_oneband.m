function [plottedLine] = PlotStatesShadedv3_oneband(Time,qS,color,transparency)
%% PLOTSTATESSHADED_V3_ONEBAND
%  Plot a nicely formatted line with the ability to tune colors and transparency
%
%  INPUTS:
%  TIME: a list of datenums
%
%  QS: a `nperiods` x 3 matrix of values, where qS(:,2) represents the mean path
%      of a variable, qS(:,1) represents a lower band, and qS(:,3) represents an
%      upper band
%
%  OPTIONAL INPUTS:
%  COLOR: An RBG color to use for the time series and bands
%
%  TRANSPARENCY: A value between 0 and 1 to set transparency of shading for bands
%
%  Created 2017-02-03 BB and ELM

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

plottedLine = plot(Time,qS(:,2),'--','Color',linecolor,'LineWidth',1.5);
hold on
fill([Time;flip(Time)],[qS(:,1); flip(qS(:,3))],.5*color,'Linestyle','none','Linestyle','none','facealpha',transparency)
plot(Time,Time*0,'k','LineWidth',.25);

xlim([Time(1), Time(end)]);
axis tight;
set(gcf,'Color','w');
set(gca,'XTick', Time(1:40:end),'XMinorTick','on');
set(gca, 'FontSize', 20);
datetick('x', 'yyyy','keeplimits', 'keepticks');
box on;
