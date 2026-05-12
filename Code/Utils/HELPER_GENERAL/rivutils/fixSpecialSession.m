function fixSpecialSession(session)

global DATA CUR FILES


if strcmp(session,'w230700')
  % 
  % Sob story:  The channel 13 was discovered off in the middle of file 009, and turned back on
  % at the beginning of 010.  Data analysis revealed that it turned off in the middle of file
  % 003.  The last three observation periods of this file are therefore bogus.
  TOTOBS = 222;
  nobs = length(DATA.multi.data)
  if nobs ~= TOTOBS
	fprintf('w230700 fixSpecialSession: Something is wrong...\n');
  else
     %	
     % Here we insert bogus 'columns' of zeros so that we can access the good data normally.
     % This means that anything that uses the multi file will need to check for such zeros.	
     %	
     start_fix_obs = 53;	
     stop_fix_obs  = 148;	

     fprintf('Inserting columns of zeros between rows %d and %d...',start_fix_obs, stop_fix_obs);
     for i=start_fix_obs:stop_fix_obs	
       d = DATA.multi.data{i};
       first_cols = d(:,1:12);
       last_cols  = d(:,13:14);
       insert     = zeros(length(d(:,1)),1);
       DATA.multi.data{i} = [first_cols insert last_cols];
     end
     fprintf('done.\n');
	
  end
elseif strcmp(session,'w210700')
  % 
  % Sob story:  Channel 13 was not turned on until file 3.

  TOTOBS = 282;
  nobs = length(DATA.multi.data)
  if nobs ~= TOTOBS
	fprintf('w210700 fixSpecialSession: Something is wrong...\n');
  else
     %	
     % Here we start off with bogus 'columns' of zeros at the end, and need to move them in a bit, 
     % to the third from the end.
     % This means that anything that uses the multi file will need to check for such zeros.	
     %	
     start_fix_obs = 1;	
     stop_fix_obs  = 101;	

     fprintf('Moving column of zeros between rows %d and %d...',start_fix_obs, stop_fix_obs);
     for i=start_fix_obs:stop_fix_obs	
       d = DATA.multi.data{i};
       first_cols = d(:,1:10);
       last_cols  = d(:,11:12);
       insert     = zeros(length(d(:,1)),1);
       DATA.multi.data{i} = [first_cols insert last_cols];
     end
     fprintf('done.\n');
	
  end

end

