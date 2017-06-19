function y = year(d,f)
%YEAR Year of date.
%   Y = YEAR(D) returns the year of a serial date number or a date string, D.
%
%	Y = YEAR(S,F) returns the year of one or more date strings S using
%   format string F. S can be a character array where each
%	row corresponds to one date string, or one dimensional cell array of
%	strings.
%
%	All of the date strings in S must have the same format F, which must be
%	composed of date format symbols according to Table 2 in DATESTR help.
%	Formats with 'Q' are not accepted.
%
%
%   For example, y = year(728647) or y = year('19-Dec-1994') returns y = 1994.
%
%   See also DATEVEC, DAY, MONTH.

%   Copyright 1995-2009 The MathWorks, Inc.

if nargin < 1
  error(message('finance:year:missingInputs'))
end

if nargin < 2
  f = '';
end

tFlag = false;   %Keep track if input was character array
if ischar(d)
  d = datenum(d,f);
  tFlag = true;
end

% Generate date vectors
if nargin < 2 || tFlag
  c = datevec(d(:));
else
  c = datevec(d(:),f);
end

y = c(:,1);             % Extract years
if ~ischar(d)
  y = reshape(y,size(d));
end

