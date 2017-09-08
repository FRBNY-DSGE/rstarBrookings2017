function bands(x,varargin)
if nargin<10
    transparency = .5;
end;
if nargin < 9
r=0.8;
g=0.8;
b=0.8;
color = [r g b];
linecolor = 'k';
else
    linecolor = color*.65;
end;

nvarargin = length(varargin);
for k = 1:nvarargin
    if k == 2
        color = [0.925    0.325    0.100]; % orange
%         color = [0.9290    0.6940    0.1250]; % gold
%         color = [0.4940    0.1840    0.5560]; % purple
%         color = [0.4660    0.6740    0.1880]; % green
%         color = [0.3010    0.7450    0.9330];
%         color = [0.6350    0.0780    0.1840]; % red
        linecolor = 0.9*color;
    end
    y = varargin{k};
    nbands = (size(y,2)-1)/2;
    for l = nbands:-1:0
        if l == 0
            plot(x,y(:,nbands+1),'--','Color',linecolor,'LineWidth',1.5);
            hold on;
        else
            (k*4/3-1/3)*0.75*color
        fill([x;flip(x)],[y(:,l); flip(y(:,end+1-l))],(1+(k-1)*1/4)*0.75*color,'LineStyle','none','FaceAlpha',transparency/sqrt(k));
        hold on;
        end
    end
end
plot(x,zeros(size(x)),'k');

end