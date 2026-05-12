
global DATA CUR FILES


% w210700
% 
% 
% Sob story:  Streamer Channel 13 was not turned on until file 3.
%
TOTOBS = 282;

for filenum = 1:3
  MUA = {};	
  file = sprintf('w210700_rivalry_00%d.adf',filenum);
  fprintf('Reconstructing abnormal file %s\n',file);
  adffile = sprintf('%s/%s',FILES.ADFPath,file)
  [p,name,e,v] = fileparts(adffile);
  loadDataFile(sprintf('%s.dgz',name));

  [nchan,nobs,sampt] = adf_getInfo(adffile);
  for obs=1:nobs
    fprintf('OBS %d\n',obs);
    for j=1:10
       j	
      [MUA.data{obs}(:,j), resampt] =  readMultiObs(adffile,j,obs);
    end
    %- here DUPLICATE last channel to avoid nan's etc in the future...	
    [MUA.data{obs}(:,11), resampt] =   readMultiObs(adffile,10,obs);
    for j=12:13
 	j
      [MUA.data{obs}(:,j), resampt] =  readMultiObs(adffile,j-1,obs);
    end
 end
 MUA.resampt = resampt;
 [p,name,e,v] = fileparts(adffile);
 muafile = sprintf('%s/%s.%s', FILES.PROC_MULTIPath, name, FILES.MULTIExt);
 save(muafile,'MUA');
end


%%%%%  ALSO!  %%%%%

for filenum = 1:3
  LFP = {};	
  file = sprintf('w210700_rivalry_00%d.eeg',filenum);
  fprintf('Reconstructing abnormal file %s\n',file);
  eegfile = sprintf('%s/%s',FILES.EEGPath,file)
  [p,name,e,v] = fileparts(eegfile);
  loadDataFile(sprintf('%s.dgz',name));
  
  [nchan,nobs,sampt] = adf_getInfo(eegfile);
  for obs=1:nobs
   fprintf('OBS %d\n',obs);
   if nchan == 12	
      fprintf('LFP electrode 5 replaced by 4\n');
      for j=1:3
        [LFP.data{obs}(:,j), resampt] =  readLFPObs(eegfile,j,obs);
      end
      [LFP.data{obs}(:,4), resampt] =  readLFPObs(eegfile,3,obs);
      for j=5:13
        [LFP.data{obs}(:,j), resampt] =  readLFPObs(eegfile,j-1,obs);
      end
   else
      fprintf('No replacement...\n');
      for j=1:nchan
        [LFP.data{obs}(:,j), resampt] =  readLFPObs(eegfile,j,obs);
      end
   end	
  end
  LFP.resampt = resampt;	

 [p,name,e,v] = fileparts(eegfile);
  lfpfile = sprintf('%s/%s.%s', FILES.PROC_LFPPath, name, FILES.LFPExt);
  save(lfpfile,'LFP');
end



