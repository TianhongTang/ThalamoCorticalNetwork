newadfdat = fixADFData(dgz,adfdata,sampt)

adf_nobs = length(adfdata);
dgz_nobs = length(dgz.e_times);

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
    adf_len = round(length(adfdata{indx})*sampt);
    difflen = dgz_len-adf_len;	
    if difflen > 5 
	next_adf_len = round(length(adfdata{indx+1})*sampt);
	if (dgz_len - adf_len - next_adf_len) < 200
          fprintf('FOUND TRUNCATED PERIOD!, DGZ OBSPER %d\n',i);
          fprintf('COMBINING MULTI OBSPERS %d and %d\n\n', indx, indx+1);
          % Here there combination of two observation periods, and a resulting difference 
          % in the increment between the dgz-related stuff and the adf-related stuff.	
          data = [adfdata{indx};adfdata{indx+1}];
          indx = indx+1;	
        end
    else
       data = adfdata{indx};
    end	
    tmpmulti{i}      = data;
  end
  newadfdata    = tmpmulti;
end	