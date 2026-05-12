function [allmua,resampt] = loadMUA(filelist,obs)

global FILES 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loadMUA  -- reads in a .mua file made previously		         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2,  obs = [];  end

allmua = {};

if isempty(obs)
  if iscell(filelist)
    len = length(filelist);
    for i=1:len
      s = muaProcessLoad(filelist{i},obs);
      if i==1
        allmua = s.MUA.data;
      else
        allmua = muaAppend(allmua,s.MUA.data);
      end
    end
  else
      
    fullpath = filelist;
    
    [p,name,e,v] = fileparts(fullpath);
    fullpath = sprintf('%s/%s.mua', FILES.PROC_MULTIPath,name);
    d = dir(fullpath);
    if ~length(d)
      adffile   = sprintf('%s/%s.adf', FILES.ADFPath, name);
      muafile 	= sprintf('%s/%s.mua', FILES.PROC_MULTIPath, name);
      s =  muaProcessLoad(fullpath,[]);
    end	
    s = load(fullpath,'-MAT');
    dgzfile  = sprintf('%s/%s.dgz',FILES.DGZPath,name);
    dgz = dg_read(dgzfile);
    allmua = muaAppend(allmua,s.MUA.data);

  end
else 
  % get only one MUA observation period only
  if iscell(filelist)
     % if combining from many files, must select the correct obs period...
     [singledgzpath,sobs] = getConcatObs(obs);	
     [p,name,e,v] = fileparts(singledgzpath);
     muapath = sprintf('%s/%s.%s', FILES.PROC_MULTIPath,name,FILES.MULTIExt);
     s =  muaProcessLoad(muapath,sobs);
   else	
     fullpath = filelist;
     [p,name,e,v] = fileparts(fullpath);
     muapath = sprintf('%s/%s.%s', FILES.PROC_MULTIPath,name,FILES.MULTIExt);
     s =  muaProcessLoad(muapath,obs);
  end
  allmua = s.MUA.data;
end
resampt = s.MUA.resampt; 

%--------------------
function s = muaProcessLoad(fullpath,obs)
global FILES

[p,name,e,v] = fileparts(fullpath);
d = dir(fullpath);
if ~length(d)
  fprintf(' loadMUA: Processing %s*\n', name);
  adffile = sprintf('%s/%s.%s', FILES.ADFPath, name, FILES.ADFExt);
  muafile = sprintf('%s/%s.%s', FILES.PROC_MULTIPath, name, FILES.MULTIExt);
  processMUAFile(adffile,muafile);
end

s = load(fullpath,'-MAT');
 if ~isempty(obs)
    if isfield(s,'MUA') & isfield(s.MUA,'data')
      s.MUA.data = s.MUA.data{obs};  % need just this observation 
      s.MUA.resampt = s.MUA.resampt;  % in sec
    else	
      s.MUA.data = s.MUL.data{obs};  % need just this observation 
      s.MUA.resampt = s.MUL.resampt;  % in sec
    end
end

%--------------------
function allmua = muaAppend(allmua,newmua)

adf_nobs = length(newmua);
all_len = length(allmua);
new_len = length(newmua);
for i=1:new_len
   allmua{all_len+i} = newmua{i};	
end
