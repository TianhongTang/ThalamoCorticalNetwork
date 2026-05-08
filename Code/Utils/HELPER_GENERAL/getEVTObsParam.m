function retlist = getEVTObsParam(dgzdata, obs, evtype, prmindx, evsubtype)
%
% Read event parameters from a dgz file
% DAL JAN-2000
%

global EVT

if (nargin < 4)
  error('usage: evtobsparam = getEVTObsParam(dgzdata,obs, evtype,prmindx,[evsubtype])')
end
if (nargin < 5) 
   	evsubtype = -1;
end	

%
%  First get the main parameter field from the cell array
params   = dgzdata.e_params;
types    = dgzdata.e_types;
subtypes = dgzdata.e_subtypes;

if (evtype == EVT.STRINGS_1) | (evtype == EVT.STRINGS_3)...
	 | (evtype == EVT.STRINGS_3)
   strflag = 1;
   retlist = {};
   tot = 0;
else
   strflag = 0;
   retlist = [];
end 
if evsubtype == -1 
	%select based only on evtype
   	ttmp = types{obs};
   	ptmp = params{obs};
	for i=1:length(ttmp)
	   type  = ttmp(i);
	   plist = ptmp{i};
	   if type == evtype
	      if ~strflag
		 if (length(plist) >= prmindx)
    	            retlist = [retlist; plist(prmindx,:)]; 
	 	 else 
		    retlist = [retlist; -1001];
		 end
	      else
                 new = plist(prmindx,:);
		 if isempty(new), break; end	
		 tot = tot+1;
  	         retlist{tot} = new;
	      end
	   end 
	end	
else
	%select based on evtype and evsubtype
   	ttmp = types{obs};
   	ptmp = params{obs};
   	stmp = subtypes{obs};
	for i=1:length(ttmp)
	   type  = ttmp(i);
	   plist = ptmp{i};
	   subtype = stmp(i);
	   if type == evtype & subtype == evsubtype
	      if ~strflag	
		 if (length(plist) >= prmindx)
    	         	retlist = [retlist; plist(prmindx,:)];
	 	 else 
			retlist = [retlist; -1001];
		 end
	      else
                 new = plist(prmindx,:);
		 if isempty(new), break; end	
		 tot = tot+1;
  	         retlist{tot} = new; 
	      end
	   end 
	end	
end

