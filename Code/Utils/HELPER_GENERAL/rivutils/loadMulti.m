function [multi,resampt] = loadMulti(filelist)

global FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loadMulti  -- reads in a .multi file made previously		         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

multi = {};

if iscell(filelist)
 len = size(filelist,2);
 for i=1:len
     fullpath = filelist{i};
     %	
     % All to make the corresponding dgz available to compensate for the
     % streamer bug.
     %
     [p,name,e,v] = fileparts(fullpath);
     fullpath = sprintf('%s/%s.multi',FILES.PROC_MULTIPath,name)  
     d = dir(fullpath);
     if ~length(d)
	[p,name,e,v] = fileparts(fullpath);
	adffile   	= sprintf('%s/%s.adf', FILES.ADFPath, name);
	multifile 	= sprintf('%s/%s.multi', FILES.PROC_MULTIPath, name);
        processMultiFile(adffile,multifile);
     end	
     dgzfile  = sprintf('%s/%s.dgz',FILES.DGZPath,name);
     dgz = dg_read(dgzfile);

     fprintf('Loading multi file %s\n',fullpath);
     s = load(fullpath,'-MAT');
     [multi,resampt] = multiAppend(multi,s.MUL.data,s.MUL.resampt,dgz);
  end
else
  fullpath = filelist;
  [p,name,e,v] = fileparts(fullpath);
  fullpath = sprintf('%s/%s', FILES.PROC_MULTIPath,filelist);
  d = dir(fullpath);
  if ~length(d)
      adffile   	= sprintf('%s/%s.adf', FILES.ADFPath, name);
      multifile 	= sprintf('%s/%s.multi', FILES.PROC_MULTIPath, name);
      processMultiFile(fullpath);
  end	
  s = load(fullpath,'-MAT');
  dgzfile  = sprintf('%s/%s.dgz',FILES.DGZPath,name);
  dgz = dg_read(dgzfile);

  [multi,resampt] = multiAppend(multi,s.MUL.data,s.MUL.resampt,dgz);
end



%--------------------

function [allmulti,resampt] = multiAppend(allmulti,newmulti,resampt,dgz)

dgz_nobs = length(dgz.e_times);
adf_nobs = length(newmulti);

%
% A bit ugly, but just in case the number of channels INCREASES during the
% recordings (BY ONE!), this will ensure that all the data gets read by
% adding a column of zeros to the end.
% Of course the fixSpecialSesssion will then have to make some more repairs.
%
if ~isempty(allmulti)
  allnchan = size(allmulti{1},2);
  newnchan   =size(newmulti{1},2);
  if allnchan < newnchan
    nobs = length(allmulti);		
    for i=1:nobs
	sz = size(allmulti{i},1);		
        allmulti{i} = [allmulti{i} zeros(sz,1)];
	fprintf('Adding Zero Period to Compensate for Unequal Channels !\n');
    end   	
  end
end

tmpmulti = [];
if adf_nobs > dgz_nobs
  % 
  % In this case, the stupid streamer divided one observation periods into two.
  % Here we simply stick the physiologic data back together.	
  %	
  fprintf('\n\nWARNING: dgz_nobs = %d, adf_nobs = %d\n',dgz_nobs, adf_nobs);
  indx = 0;
  for i=1:dgz_nobs
    indx = indx+1;	
    dgz_len = round(max(dgz.e_times{i}));
    adf_len = round(length(newmulti{indx})*resampt);
    difflen = dgz_len-adf_len;	
    if difflen > 5 
	next_adf_len = round(length(newmulti{indx+1})*resampt);
	if (dgz_len - adf_len - next_adf_len) < 200
          fprintf('FOUND TRUNCATED PERIOD!, DGZ OBSPER %d\n',i);
          fprintf('COMBINING MULTI OBSPERS %d and %d\n\n', indx, indx+1);
          % Here there combination of two observation periods, and a resulting difference 
          % in the increment between the dgz-related stuff and the adf-related stuff.	
          data = [newmulti{indx};newmulti{indx+1}];
          indx = indx+1;	
        end
    else
       data = newmulti{indx};
    end	
    tmpmulti{i}      = data;
  end
  newmulti    = tmpmulti;
elseif adf_nobs < dgz_nobs
  fprintf('WARNING! MORE DGZ OBSPERS THAN ADF OBSPERS!\n');
  % 
  % This can only be for stupid reasons--e.g. that the channel was turned off in the 
  % middle of the file. For now, just add some null obspers and deal with it in the
  % fix file.
  %	
  for i=1:dgz_nobs
    dgz_len = round(max(dgz.e_times{i}));
    if i <= adf_nobs	
      adf_len = round(length(newmulti{i})*resampt);
      next_dgz_len = round(max(dgz.e_times{i+1}));
      difflen = dgz_len-adf_len;	 
      %    fprintf('\n----------------------------------\n');
      %    fprintf('N DGZ OBSP = %d, N ADF OBSP = %d\n',dgz_len, adf_len);
      if abs(difflen) > 5 
	% Here set to null any potentially bogus observation period data sets.  We will
	% check for such null periods later.
        data = [];
      else	
        data = newmulti{i};
      end
    else 
      data = [];      	 	     
      fprintf('ADDING NULL DATA PERIOD!\n');
    end
    tmpmulti{i}      = data;
  end
  newmulti    = tmpmulti;
end

all_len = length(allmulti);
new_len = length(newmulti);

for i=1:length(newmulti)
   allmulti{all_len+i} = newmulti{i};	
end

