function bands(x,varargin)

% vars = {'plotLine', 'transparency', 'color', 'linecolor'};
% defaults = {1, 0.5, [0.8 0.8 0,8], 'k'};
% varargparse(varargin, vars, defaults);


if nargin< 10
  transparency = .5;
else
  transparency = varargin{9};
end;
if nargin < 9
  r=0.8;
  g=0.8;
  b=0.8;
  color = [r g b];
  linecolor = 'k';
else
  color = varargin{8};
  linecolor = color*.65;
end;

nvarargin = length(varargin);
for k = 1:nvarargin
    if k == 2
        color = [0.925    0.325    0.100]; % orange
        linecolor = 0.9*color;
    end
    y = varargin{k};

    % If number of columns is odd -> it contains mean or median line
    if mod(size(y,2),2) == 1
      containsMedianLine = 1;
      nbands = (size(y,2)-1)/2;
    else
      containsMedianLine = 0;
      nbands = (size(y,2))/2;
    end
    for l = nbands:-1:0
        if l == 0
	  if containsMedianLine
            plot(x,y(:,nbands+1),'--','Color',linecolor,'LineWidth',1.5);
            hold on;
	  end
        else
            (k*4/3-1/3)*0.75*color;
	  fill([x;flip(x)],[y(:,l); flip(y(:,end+1-l))],(1+(k-1)*1/4)*0.75*color,'LineStyle','none','FaceAlpha',transparency/sqrt(k));
	  hold on;
        end
    end
end
plot(x,zeros(size(x)),'k');

set(gcf,'Color','w')
set(gca,'XTick', x(1:40:end),'XMinorTick','on')
axis tight; box on;
datetick('x', 'yyyy', 'keeplimits', 'keepticks')

end
