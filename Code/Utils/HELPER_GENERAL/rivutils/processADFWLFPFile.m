function LFP = processADFWLFPFile(adfwfile,lfpfile)

  [nchan,nobs,sampt] = adfw_info(adfwfile);
  wb = waitbar(0, sprintf('Processing LFP Data'));
  for obs=1:nobs
    waitbar((obs/nobs),wb);
    for j=1:nchan
      [LFP.data{obs}(:,j), resampt] =  readADFWLFPObs(adfwfile,j,obs);
    end
  end
  LFP.resampt = resampt;	
  close(wb);
  
  %
  % Here just save the thing for future use
  save(lfpfile,'LFP');

