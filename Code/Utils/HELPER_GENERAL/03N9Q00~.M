function retlist = getEVTParams(dgzdata, evtype, prmindx, evsubtype)
%
% Read event parameters from a dgz file
% DAL JAN-2000
% 16-Aug-2000  YM modified from DAL version, using find()
global EVT

if (nargin < 3)
  error('usage: evtparams = getEVTParams(dgzdata,evtype,prmindx,[evsubtype])')
end
if (nargin < 4), evsubtype = -1;  end	

OLDSTRINGS = 132;
if (evtype == EVT.STRINGS_1) | (evtype == EVT.STRINGS_2) ...
	    | (evtype == EVT.STRINGS_3) | (evtype == OLDSTRINGS) 
  strflag = 1;
else
  strflag = 0;
end


% select events
tmpsel = {};
if evsubtype == -1 
  %select based only on evtype
  for j=1:length(dgzdata.e_types)
    tmpsel{j} = find(dgzdata.e_types{j} == evtype);
  end
else
  %select based on evtype and evsubtype
  for j=1:length(dgzdata.e_types)
    tmpsel{j} = find(dgzdata.e_types{j} == evtype &...
	       dgzdata.e_subtypes{j} == evsubtype);
  end
end

% get parameters
retlist = {};
if ~strflag
  for j=1:length(tmpsel)
    if length(tmpsel{j}) ~= 0
      rettmp = [];
      plistS = dgzdata.e_params{j}(tmpsel{j});
      for i=1:length(plistS)
	plist = plistS{i};
	if length(plist) >= prmindx
	  rettmp = [rettmp; plist(prmindx)];
	else
	  rettmp = [rettmp; -1001];
	end
      end
      retlist{j} = rettmp;
    else
      retlist{j} = [];
    end
  end
else
  for j=1:length(tmpsel)
    if length(tmpsel{j}) ~= 0
      rettmp = {};
      plistS = dgzdata.e_params{j}(tmpsel{j});
      for i=1:length(plistS)
	plist = plistS{i};
	if size(plist,1) >= prmindx
	  str = deblank(plist(prmindx,:));
	else
	  str = '';
	end
	rettmp{i} = str;
      end
      retlist{j} = rettmp;
    else
      retlist{j} = {};
    end
  end
end
