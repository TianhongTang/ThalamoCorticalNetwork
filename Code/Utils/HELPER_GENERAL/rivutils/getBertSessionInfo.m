function SI = getSessionInfo(file)
SI.angle 	= [];
SI.rf 		= [];
SI.spontfile	= [];
fr = getFileRoot(file);
frr = fr(1:10);
switch(frr)
  case {'bert_18_04_2_'}
    SI.electrodes = [ 1  2 3 5 6 7 8 9 10 11 15];
    SI.areas      = [-1 -1 4 4 4 4 4 4  4  2  2]; 
    SI.spontfile  = 2;  % REM-sleep?
  case {'bert_19_04_02'}
    SI.electrodes = [ 1  2 3 5 6 7 8 9 10 11 15];
    SI.areas      = [-1 -1 4 4 4 4 4 4  4  2  2]; 
    SI.spontfile  = 2;  % REM-sleep?
end


