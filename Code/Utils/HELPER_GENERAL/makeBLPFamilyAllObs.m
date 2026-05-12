function BLP = makeBLPFamilyAllObs(alldat,dumpfile,Fs)

  nobs  = size(alldat);
  nchan = size(alldat{1},2);
  
  decval = 20;
  decval1 = 4;
  decval2 = decval/decval1;

  frange = getFrequencyRanges;

  for obs = 1:nobs
	fprintf('OBS %d:',obs);
	for chan=1:nchan
	  fprintf(' %d',chan);
	  dec1trace = decimate(alldat{obs}(:,chan),decval1);
	  dec1Fs = Fs/decval1;

	  [dat,Fs_BLP] =  makeBLPTrial(dec1trace,dec1Fs,decval2,frange);
	  if chan == 1
		BLP.dat{obs} = zeros(size(dat,1),size(dat,2),nchan);
	  end
	  BLP.dat{obs}(:,:,chan) = dat;
	end
	fprintf('...done\n');;
  end
  BLP.Fs = Fs_BLP;
  BLP.frange = frange;

  save(dumpfile,'BLP');
  
  
function [blp_dat,Fs_BLP] = makeBLPTrial(dat,Fs,decval,frange)
  
  
  nbw    = size(frange,2);
  nyq = Fs/2;		    	                    % Hz

  for bw = 1:nbw
    fr            = frange(:,bw);
    hpc 		  = fr(1);	    	            % Hz
    lpc 		  = fr(2);	    	            % Hz
	if lpc > nyq, continue; end
    hWn  		  = hpc/nyq;		   	        % between 0.0 and 1.0
    lWn  		  = lpc/nyq;   	   	            % between 0.0 and 1.0
    [lb,la]       = cheby1(2,0.8,[hWn lWn]);	% chebyshev bandpass
    fdat          = filtfilt(lb,la,dat);
    blp_dat(:,bw) = decimate(abs(fdat),decval); % bandpass low
  end
  blp_dat  = blp_dat;
  Fs_BLP   = Fs/decval;
    
  
  
  
function franges = getFrequencyRanges()
  
  franges = [9 14;...
	     15 30;...
	     30 50;...
	     50 100]';
  