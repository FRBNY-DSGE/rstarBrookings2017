function [out] = varargparse(args, params, defaults)
%VARARGPARSE Parse varargin given expected parameters and defaults
%
%   VARARGPARSE(ARGS, PARAMS, DEFAULTS) parses ARGS given expected parameters
%   PARAMS and default values DEFAULTS. Param-value pairs that appear in ARGS
%   but not in PARAMS cause errors. PARAMS and DEFAULTS are both cell arrays
%   where a parameter in PARAMS has a default value in the same position in
%   DEFAULTS.
%
%   Due to the way MATLAB resolves variable names, if any of the names in PARAMS
%   is a built in function, you may run into problems. Resolve this by adding
%   initializations before you call VARARGPARSE.
%
%   There are many, many ways to process optional arguments. The built-in
%   INPUTPARSER works well but is quite verbose. See more discussion on
%   http://stackoverflow.com/q/2775263/2514228.
%
%   Example
%       % Define function that takes varargs
%       function [] = foo(varargin)
%       params = {'baz','qux'};
%       defaults = {1, 3};
%       varargparse(varargin, params, defaults);
%       disp(baz)
%       disp(qux)
%       end
%
%       % Call function with param-value pairs.
%       foo('bar', 2);
%
%   See also INPUTPARSER, VALIDATESTRING, VALIDATEATTRIBUTES

% Created 2016-04-13 MJS
expected = params;

% Ensure we have a matching number of params and values.
if mod(length(args),2) ~= 0
  error('Invalid param-value arguments.');
end

% Extract provided params (odd elements) and values (even elements)
params   = args(1:2:end);
values   = args(2:2:end);

% For each expected param, either pop it from the provided params or set it to
% its default value.
for i=1:length(expected);
  key = expected{i};
  if ismember(key, params)
    j = find(strcmpi(key, params));
    assignin('caller', key, values{j});
    params(j) = [];
    values(j) = [];
  else
    assignin('caller', key, defaults{i});
  end
end

% Any additional param-value pairs are invalid.
if length(params) > 0
  error(['Invalid param(s): ', params{:}]);
end

end % of varargparse
