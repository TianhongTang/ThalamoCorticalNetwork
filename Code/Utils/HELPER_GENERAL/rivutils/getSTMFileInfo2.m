function stm = getSTMFileInfo2(stmfile)
if nargin < 1
  stm = getSTMFileInfo;
else
  stm = getSTMFileInfo(stmfile);
end
