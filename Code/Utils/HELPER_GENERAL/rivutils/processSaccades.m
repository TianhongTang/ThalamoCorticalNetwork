function processSaccades(file,dir)
%
%  Process saccades/fixation periods.  
%

global FILES  
if ~nargin  
  loadDataFile
else
file
dir
  loadDataFile(file,dir)
end
emFileProcess(FILES.FILERoot,FILES.DGZPath,FILES.PROC_SACPath);


































