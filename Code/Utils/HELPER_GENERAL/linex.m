function [h] = linex(val,varargin)
% [h] = linex(val,varargin)
%
% plots a line that extends for a long way along the y axis
%  if given a vector, plots one line for each entry
%
% 'for a long way' means 'until 10^12'

h = plot_line_xint(val,varargin{:});
