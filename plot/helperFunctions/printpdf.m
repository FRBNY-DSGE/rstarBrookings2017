function printpdf(h,outfilename, varargin)

vars = {'square', 'fontsize'};
defaults = {0, 20};
varargparse(varargin, vars, defaults);

axs = findobj(h, 'type', 'axes');
for ax = axs
  set(ax,'TickLabelInterpreter', 'latex')
end

if square
  pbaspect([1,1,1]);
else
  pbaspect([( 1 + sqrt(5) ) / 2, 1, 1])
end

set(h, 'PaperUnits', 'inches');
set(h, 'Units', 'inches');
pos=get(h,'Position');

set(gca, 'fontsize', fontsize);

if square
  set(h,'PaperSize',[pos(3) pos(3)]);
  set(h,'PaperPositionMode', 'manual');
  set(h, 'PaperPosition',[0 0 pos(3) pos(3)]);
else
  set(h,'PaperSize',[pos(3) pos(4)]);
  set(h,'PaperPositionMode', 'manual');
  set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
end

print('-dpdf', outfilename);
print('-depsc', outfilename);

end