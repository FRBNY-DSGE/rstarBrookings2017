function [barPos barNeg l1] = plotStackedShockdec(time, shockdec, mean, ...
    shockCategories, barColors, barLabels, varargin)
%% PLOTSTACKEDSHOCKDEC
%  Plot a stacked shock decomposition, with shocks grouped into categories.
%
%  INPUTS:
%  TIME:     a list of datenums for the dates the shockdec covers
%
%  SHOCKDEC: a `nperiods` x `nshocks` database representing the contributions of
%            each shock to the deviation of a variable from steady-state
%
%  MEAN:     the detrended mean of a time series represending the path of a
%            variable
%
%  SHOCKCATEGORIES: a cell array of integer arrays. Each integer represents
%                   the index of a shock in SHOCKDEC to be included in that
%                   category. Each bar in the resulting plot represents the
%                   cumulated impact of the shocks in that category.
%
%  BARCOLORS: A cell array of RGB values or MATLAB color names, in the same
%             order as SHOCKCATEGORIES
%
%  BARLABELS: A cell array of labels for the bars
%
%  OPTIONAL (NAMED) INPUTS
%  FIG: a figure. Usage:
%    f = figure()
%    [barPos, barNeg, l1] = plotStackedShockdec(...'fig', f)
%
%  Created 2017-02-14 ELM

vars = {'fig'};
defaults = {'0'};
varargparse(varargin, vars, defaults);

if fig == 0
  fig = figure();
end

% Split the positive and negative parts of the shockdec
shockdecPos = shockdec;
shockdecNeg = shockdec;
shockdecPos(shockdecPos < 0) = 0;
shockdecNeg(shockdecNeg > 0) = 0;

% Plot 2 separate graphs, one positive and one negative

hold on;
barPos = bar(time, shockdecPos,'stacked','edgecolor', 'none');
barNeg = bar(time, shockdecNeg,'stacked','edgecolor', 'none');

% Plot the mean of the variable and the x-axis
l1 = plot(time, mean, 'k', 'LineWidth', 2.0);
plot(time,time*0,'k','LineWidth',.25)

% Format
axis tight;
set(gca, 'FontSize', 20);
tick10 = time(1:40:end);
set(gca, 'XTick', tick10);
datetick('x', 'yyyy','keeplimits', 'keepticks');

% Set bar colors
for iList = 1:length(barColors)
  if isnumeric(barColors{iList})
    set(barPos(iList),'FaceColor',barColors{iList});
    set(barNeg(iList),'FaceColor',barColors{iList});
  else
    set(barPos(iList),'FaceColor',rgb(barColors{iList}));
    set(barNeg(iList),'FaceColor',rgb(barColors{iList}));
  end
end

box on;
hold off;

end