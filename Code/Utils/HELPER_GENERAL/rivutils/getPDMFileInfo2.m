function pdm = getPDMFileInfo2(pdmfile)
if nargin < 1
  pdm = getPDMFileInfo;
else
  pdm = getPDMFileInfo(pdmfile);
end
