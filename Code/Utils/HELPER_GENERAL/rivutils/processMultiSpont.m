function SP = processMultiSpont()

global FILES

adffile 	= sprintf('%s/%s.adf',FILES.SPONTPath,FILES.SPONTFile);
multifile 	= sprintf('%s/%s.multi',FILES.PROC_MULTIPath,FILES.SPONTFile);
MUL = checkFile(multifile,'MUL');
if ~isempty(MUL)
  fprintf('Successfully read processed multi file %s\n', multifile);
else
  MUL = processMultiFile(adffile,multifile);
end;

SP.data = MUL.data;
SP.sampt = MUL.resampt;
