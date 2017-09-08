function printpdf(h,outfilename)

set(h, 'PaperUnits', 'inches');
set(h, 'Units', 'inches');
pos=get(h,'Position');
set(h,'PaperSize',[pos(3) pos(4)]);
set(h,'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf', outfilename);

end