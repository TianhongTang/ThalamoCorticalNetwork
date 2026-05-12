function [singledgzpath, singleobs] = getConcatObs(obs)
% VERSION : 1.00  08-Apr-01  DAL
%           1.01  10-Apr-01  YM,  checks datafiletype

global FILES

if (strcmp(FILES.datafiletype,'CONCATFILE'))
  filelist = FILES.curfiles;
  len = length(filelist);
  cumsum = 0;
  lastcumsum = 0;
  for i=1:len
    path = sprintf('%s/%s.dgz',FILES.DGZPath,filelist{i});
    tmp  = dg_read(path);
    cumsum = cumsum+length(tmp.e_types);
    if obs > lastcumsum & obs <= cumsum
      singleobs  = obs-lastcumsum;
      singledgzpath = path;
      break;	
    end
    lastcumsum = cumsum;	
  end
else
  singleobs     = obs;
  singledgzpath = sprintf('%s/%s',FILES.DGZPath,FILES.DGZFile);
end

