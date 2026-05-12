function SI = getSessionInfo(file)
SI.angle 	= [];
SI.rf 		= [];
SI.spontfile	= [];
fr = getFileRoot(file);
frr = fr(1:6);

switch(frr)
  case {'030802'}
    SI.electrodes = [ 3 4 5 6 7 9 11];
    SI.areas      = [ 0 0 0 0 0 0 0]
  case {'040802'}
    SI.electrodes = [ 4 5 6 7 9 11];
    SI.areas      = [ 0 0 0 0 0 0]
  case {'050802'}
    SI.electrodes = [ 5 6 7];
    SI.areas      = [ 0 0 0]
  case {'060802'}
    SI.electrodes = [ 4 5 6 7 8];
    SI.areas      = [ 0 0 0 0 0]
  case {'080802'}
    SI.electrodes = [ 4 6 7];
    SI.areas      = [ 0 0 0]
  case {'090802'}
    SI.electrodes = [ 4 6 7 12 15];
    SI.areas      = [ 0 0 0 0 0]
  case {'100802'}
    SI.electrodes = [ 4 5 6 11 15];
    SI.areas      = [ 0 0 0 0 0]
  case {'110802'}
    SI.electrodes = [ 4 5 7 12 15];
    SI.areas      = [ 0 0 0 0 0]
  case {'120802'}
    SI.electrodes = [ 4 7 12 15];
    SI.areas      = [ 4 4 5 5]
  case {'130802'}
    SI.electrodes = [ 4 5 6 7 9 12 15];
    SI.areas      = [ 0 0 0 0 0 0 0];
  case {'140802'}
    SI.electrodes = [ 4 5 6 7 9 12 15];
    SI.areas      = [ 0 0 0 0 0 0 0];
  case {'150802'}
    SI.electrodes = [ 4 6 9 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0];
  case {'160802'}
    SI.electrodes = [ 4 5 7 8 9 11 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'170802'}
    SI.electrodes = [ 4 5 7 8 9 11 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'190802'}
    SI.electrodes = [ 4 5 7 8 9 11 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'200802'}
    SI.electrodes = [ 4 5 7 8 9 11 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'020902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'030902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'040902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'050902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'060902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];  
  case {'070902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'080902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'090902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'100902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'110902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];  
  case {'120902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];  
  case {'130902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'140902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'150902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'160902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];
  case {'170902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0]; 
  case {'180902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0]; 
  case {'200902'}
    SI.electrodes = [ 4 10 12 15];
    SI.areas      = [ 0 0 0 0]; 
  case {'210902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0]; 
  case {'220902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0]; 
  case {'230902'}
    SI.electrodes = [ 2 4 6 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0]; 
  case {'240902'}
    SI.electrodes = [ 4 5 6 7 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0 0]; 
  case {'wally_'}
    SI.electrodes = [ 4 5 6 7 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0 0]; 
  case {'250902'}
    SI.electrodes = [ 4 5 6 7 8 10 12 14 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0 0]; 
  case {'290902'}
    SI.electrodes = [ 4 5 7 10];
    SI.areas      = [ 0 0 0 0]; 
  case {'300902'}
    SI.electrodes = [ 4 5 7 10];
    SI.areas      = [ 0 0 0 0]; 
  case {'011002'}
    SI.electrodes = [ 4 5 7 10];
    SI.areas      = [ 0 0 0 0];  
  case {'021002'}
    SI.electrodes = [ 5 6 8 10 12 14];
    SI.areas      = [ 0 0 0 0 0 0];  
  case {'031002'}
    SI.electrodes = [ 4 5 6 10];
    SI.areas      = [ 0 0 0 0];  
  case {'041002'}
    SI.electrodes = [ 4 5 6 10];
    SI.areas      = [ 0 0 0 0];  
  case {'051002'}
    SI.electrodes = [ 6 10 12 15];
    SI.areas      = [ 0 0 0 0];  
  case {'071002'}
    SI.electrodes = [ 4 6 7 8 9 10 12 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];  
  case {'091002'}
    SI.electrodes = [ 4 6 7 8 9 10 12 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];  
  case {'111002'}
    SI.electrodes = [ 4 6 7 8 9 10 12 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];  
  case {'131002'}
    SI.electrodes = [ 4 6 7 8 9 10 12 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0];  
  case {'151002'}
    SI.electrodes = [ 4 6 7 8 9 10 12 13 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0 0];  
  case {'231002'}
    SI.electrodes = [ 4 6 7 8 9 10 12 13 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0 0];  
  case {'251002'}
    SI.electrodes = [ 4 6 7 8 9 10 12 13 15];
    SI.areas      = [ 0 0 0 0 0 0 0 0 0];     
  case {'130203'}
    SI.electrodes = [ 2 4 5 6 8 9 10 11 15 ];
    SI.areas      = [ 0 0 0 0 0 0 0 0 0 ];
end


