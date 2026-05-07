
global DATA CUR FILES


% w230700
% 
% Sob story:  The channel 13 was discovered off in the middle of file 009, and turned back on
% at the beginning of 010.  Data analysis revealed that it turned off in the middle of file
% 003.  The last three observation periods of this file are therefore bogus.

TOTOBS = 222;
start_fix_obs = 53;	
stop_fix_obs  = 148;	

%---------------------------------
% first prepare the first file...-
%---------------------------------

MUA = {};	
setUserFilePaths
file = 'w230700_rivalry_003.adf';
adffile = sprintf('%s/%s',FILES.ADFPath,file)
[p,name,e,v] = fileparts(adffile);
loadDataFile(sprintf('%s.dgz',name));
[nchan,nadfobs,sampt] = adf_getInfo(adffile);
nchan = 15;
fprintf('Reconstructing normal part of 003\n');
for obs=1:19
   fprintf('OBS %d\n',obs);
   for j=1:nchan
     [MUA.data{obs}(:,j), resampt] =  readMultiObs(adffile,j,obs);
   end
end

%
% Insert lost observation period...by duplicating subsequent one...
%
for j=1:nchan
  [MUA.data{20}(:,j), resampt] =  readMultiObs(adffile,j,21);
end

for obs=21:22
   fprintf('OBS %d\n',obs);
   for j=1:nchan
     [MUA.data{obs}(:,j), resampt] =  readMultiObs(adffile,j,obs-1);
   end
end

%for obs=nobs-2:nobs
%   fprintf('OBS %d\n',obs);
%   for j=1:12
%     [MUA.data{obs}(:,j), resampt] =  readMultiObs(adffile,j,obs);
%   end
%   %- here DUPLICATE last channel to avoid nan's etc in the analysis...	
%   [MUA.data{obs}(:,13), resampt] =   readMultiObs(adffile,12,obs);
%   for j=14:15
%	j
%     [MUA.data{obs}(:,j), resampt] =  readMultiObs(adffile,j,obs);
%   end
%end

MUA.resampt = resampt;
[p,name,e,v] = fileparts(adffile);
muafile = sprintf('%s/%s.%s', FILES.PROC_MULTIPath, name, FILES.MULTIExt);
save(muafile,'MUA');


for filenum = 4:9
  MUA = {};	
  file = sprintf('w230700_rivalry_00%d.adf',filenum);
  fprintf('Reconstructing abnormal file %s\n',file);
  adffile = sprintf('%s/%s',FILES.ADFPath,file)
  [p,name,e,v] = fileparts(adffile);
  loadDataFile(sprintf('%s.dgz',name));

  [nchan,nobs,sampt] = adf_getInfo(adffile);
  for obs=1:nobs
    fprintf('OBS %d\n',obs);
    for j=1:12
       j	
      [MUA.data{obs}(:,j), resampt] =  readMultiObs(adffile,j,obs);
    end
    %- here DUPLICATE last channel to avoid nan's etc in the future...	
    [MUA.data{obs}(:,13), resampt] =   readMultiObs(adffile,12,obs);
    for j=14:15
 	j
      [MUA.data{obs}(:,j), resampt] =  readMultiObs(adffile,j-1,obs);
    end
 end
 MUA.resampt = resampt;
 [p,name,e,v] = fileparts(adffile);
 muafile = sprintf('%s/%s.%s', FILES.PROC_MULTIPath, name, FILES.MULTIExt);
 save(muafile,'MUA');
end



