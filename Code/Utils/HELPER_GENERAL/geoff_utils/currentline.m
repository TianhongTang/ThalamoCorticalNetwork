function linenum = currentline
% CURRENTLINE   Current line number in an M-file.
%    linenum = CURRENTLINE returns the line number on which the statement
%    appears in an M-file. If called from the command window, not in an
%    M-file, linenum is set to zero. This is intended for debugging without
%    using the MATLAB debugger interface, which is particularly useful for
%    "heisenbug" situations when using the debugger changes the nature of
%    the bug.
%
%    GKA July 2007

ST=dbstack;
if numel(ST)>1
    linenum = ST(2).line;
else
    linenum = 0;
end