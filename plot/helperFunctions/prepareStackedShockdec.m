function [barPos, barNeg, l1] =  prepareStackedShockdec(filename, trendFilename, shockdecStart, ...
	shockdecEnd, shockcatGroups, shockcatNames, shockdecColors, varargin)
%% PREPARESTACKEDSHOCKDEC
%  Read in shockdec, trend, and deterministic trend, do some preparations, and
%  plot stacked shockdec.
%
%  Specifically, the trend, deterministic trend, and shock decompositions are
%  all saved separately but are needed to produce a stacked shockdec. This
%  function combines them together to produce either a shockdec of a detrended
%  variable or the variable with the trend included. It returns the positive
%  barplot, the negative barplot, and the detrended mean line plot objects.
%
%  INPUTS
%  FILENAME: Name of file containing shockdec (contributions of each shock to a
%	variable's path)
%
%  TRENDFILENAME: Name of file containing trend (estimated constant)
%
%  SHOCKDECSTART: Start date for shockdec (as a datenum)
%
%  SHOCKDECEND: End date for shockdec (as a datenum)
%
%  SHOCKCATNAMES: Cell array of names for each category of shocks
%
%  SHOCKDECCOLORS: Cell array of rgb values or MATLAB color names for each
% 	           variable
%
%  OPTIONAL (NAMED) INPUTS
%  FIG: a figure. Usage:
%    f = figure()
%    [barPos, barNeg, l1] = plotStackedShockdec(...'fig', f)

vars = {'fig'};
defaults = {'0'};
varargparse(varargin, vars, defaults);

if fig == 0
  fig = figure();
end

% Read in shockdec and sum up all the shock categories
shockdecTable = readtable(filename);
shockcats = {};
for i = 1:length(shockcatGroups)
  inds = cellfun(@(x) find(strcmp(shockdecTable.Properties.VariableNames, x)), ...
      shockcatGroups{i}, 'UniformOutput', false);
  shockcats{i} = cell2mat(inds);
end
time          = datenum(shockdecTable.date(:));
startInd      = find(time == shockdecStart);
endInd        = find(time == shockdecEnd);
time          = time(startInd:endInd);
shockdec      = zeros(size(shockdecTable(startInd:endInd,1),1), length(shockcats));
for i = 1:length(shockcats)
  shockdec(:, i) = sum(table2array(shockdecTable(startInd:endInd, shockcats{i})), 2);
end

% Read in the trend, extracting the correct periods
trend0 = 0;
mean = shockdecTable.detrendedMean(startInd:endInd);
trendTable    = readtable(trendFilename);
trendTime     = datenum(trendTable.date(:));
trendStartInd = find(trendTime == shockdecStart);
trendEndInd   = find(trendTime == shockdecEnd);
trend         = trendTable.mean(startInd:endInd);
trend0        = trend(1);

% separate out the deterministic trend
meanDettrend  = shockdecTable.dettrend(startInd:endInd);

% Add dettrend to the shockdec in the "Other" category. If "Other" category
% doesn't exist, make a category for it
shockcatOther = find(strcmp(shockcatNames, 'Other'));
if isempty(shockcatOther)
  shockdec(:, end+1) = meanDettrend;
  % shockdecColors{end+1} = 'LightGray';
  shockdecColors{end+1} = shockdecColors{end};
  shockcatNames{end+1}  = 'dt';
  shockcatOther = size(shockdec,2);
else
  shockdec(:, shockcatOther) = shockdec(:, shockcatOther) + ...
      meanDettrend;
end

% plot
[barPos, barNeg, l1] = plotStackedShockdec(time, shockdec, mean, shockcats, ...
    shockdecColors, shockcatNames, 'fig', fig);
end