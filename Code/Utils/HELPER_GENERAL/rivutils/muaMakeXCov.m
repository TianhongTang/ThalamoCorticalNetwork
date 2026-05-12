function sxcv = muaMakeXCov(file,obs,spont)

global FILES CUR DATA

ffile 		= strrep(file,'_','\_');

lags		= 128;
win 		= 256;
method 		= 'xcov';
shuffle		= 0;
maxpts		= 1000;
retval 		= 1;

if exist('spont')==1 & ~isempty(spont) & spont
  loadfile = FILES.SPONTFile;
  rootfile = getFileRoot(FILES.SPONTFile);
else 
  loadfile = FILES.MULTIFile;
  rootfile = getFileRoot(file);
end

% load XCV data
muarivcovfile = sprintf('%s_OBS%d_MUA_WN%d_LG%d.%s',...
			rootfile,obs,win,lags,method);
sxcv = checkFile(sprintf('%s/%s/%s',FILES.PROC_XCVPath,FILES.FILERoot,...
	muarivcovfile));

if isempty(sxcv)
  % load MUA data
  DATA.multi = [];
  [DATA.multi,DATA.multi_resampt] = loadMUA(loadfile,obs);
  [DATA.multi,A,E] = sortByArea(file,DATA.multi);
  sampt = DATA.multi_resampt;	
  if sampt > 0.1
     sampt = sampt/1000;
  end
  sxcv  = processXCov(DATA.multi,sampt,method,win,lags,...
		      shuffle,muarivcovfile,A,E,maxpts);  
  DATA.multi = [];
end
