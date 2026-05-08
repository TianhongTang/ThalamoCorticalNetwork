function retlist = getParams(dgzdata, evtype, prmindx, evsubtype)

if (nargin < 3)
	error('usage: getParams(dgzdata,evtype,prmindx,[evsubtype])')
end
if (nargin < 4) 
   	evsubtype = -1;
end	

%
%  First get the main parameter field from the cell array
params   = getfield(dgzdata, 'e_params');
types    = getfield(dgzdata, 'e_types');
subtypes = getfield(dgzdata, 'e_subtypes');

retlist = [];
if evsubtype == -1 
	%select based only on evtype
	for j=1:length(types)
	   	ttmp = types{j};
	   	ptmp = params{j};
		for i=1:length(ttmp)
		   type  = ttmp(i);
		   plist = ptmp{i};
		   if type == evtype
			retlist = [retlist;plist(prmindx+1)];  %note the index conversion to matlab
		   end 
		end	
	end
else
	%select based on evtype and evsubtype
	%select based only on evtype
	for j=1:length(types)
	   	ttmp = types{j};
	   	ptmp = params{j};
	   	stmp = subtypes{j};
		for i=1:length(ttmp)
		   type  = ttmp(i);
		   plist = ptmp{i};
		   subtype = stmp(i);
		   if type == evtype & subtype == evsubtype
			retlist = [retlist;plist(prmindx+1)];  %note the index conversion to matlab
		   end 
		end	
	end
end

