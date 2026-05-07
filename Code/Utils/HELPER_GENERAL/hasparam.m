function [result] = hasparam(cellparams,paramname)
% [result] = hasparam(cellparams,paramname)
%  return true if the cell array 'cellparams'
%   contains the string 'paramname', followed
%   by at least one other argument
%
%  is case-insensitive

for c = 1:(length(cellparams)-1)
    if isstr(cellparams{c}) && strcmp(lower(cellparams{c}),lower(paramname))
        result = true;
        return;
    end;
end;

result = false;
