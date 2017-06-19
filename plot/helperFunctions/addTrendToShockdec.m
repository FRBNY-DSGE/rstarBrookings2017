function [ax] = addTrendToShockdec(ax, trend, accuracy, varargin)

  vars = {'ymin', 'ymax', 'forceLower', 'forceUpper'};
  defaults = {0, 5, 0, 0};
  varargparse(varargin, vars, defaults);

  % value I want to center the shockdec around (trend rounded to nearest 0.5)
  center = round(trend/accuracy)*accuracy;

  % Step 1: Change labels to only be at 0.5 increments
  L = get(ax, 'YLim');

  % round lower axis down to nearest multiple of 0.5, or
  % the negative of the trend (so our final plot starts at 0)
  % lowerTmp = min(L(1), -trend);
  lower = min(floor(L(1) * (1/accuracy)) *accuracy, ymin -center);
  if forceLower
    lower = ymin - center;
  end

  % round upper limit up to nearest multiple of 0.5
  upper = max(ceil(L(2) * (1/accuracy)) * accuracy, ymax - center);
  if forceUpper
    upper = ymax - center;
  end

  % change YLim to go from lower to upper
  set(ax, 'YLim', [lower upper]);

  % set tempTicks to be at 0.5 incremebts
  tempTicks = lower:accuracy:upper;
  set(ax, 'YTick', tempTicks);

  % Step 2: Change the levels of the labels to center the plot around the trend
  % value rounded to the nearest 0.5

  % reassign YTick and YTickLabel.
  % - YTickLabels are reassigned to reasonable numbers centered at roughly the
  %   trend value.
  % - YTicks are reassigned to numbers that maintain the appropriate distance from
  %   the YTickLabels.

  ticks  = get(ax, 'YTick');
  labels = get(ax, 'YTickLabel');
  newTicks  = ticks + center - trend;

  newLabels = cell(size(ticks));
  for i = 1:length(ticks)
    newLabels(i) = num2cell(ticks(i) + center);
  end

  set(ax, 'YTickLabel', newLabels);

end