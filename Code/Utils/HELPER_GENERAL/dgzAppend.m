function dgzdata = dgzAppend(dgzdata,dgztmp)

if ~length(dgzdata)
   dgzdata.e_times = [];	
   dgzdata.e_types = [];	
   dgzdata.e_subtypes = [];	
   dgzdata.e_params = [];	
   dgzdata.ems = [];	
%   dgzdata.spk = {};	
   dgzdata.spk_types = [];	
   dgzdata.spk_channels = [];	
   dgzdata.spk_times = [];	
   dgzdata.obs_times = [];	
end

old_obsnum = length(dgzdata.e_times);
cur_obsnum = length(dgztmp.e_times);
new_obsnum = old_obsnum+cur_obsnum;

dgzdata.e_times      = [dgzdata.e_times;dgztmp.e_times];
dgzdata.e_types      = [dgzdata.e_types;dgztmp.e_types];
dgzdata.e_subtypes   = [dgzdata.e_subtypes;dgztmp.e_subtypes];
dgzdata.e_params     = [dgzdata.e_params;dgztmp.e_params];
dgzdata.ems          = [dgzdata.ems;dgztmp.ems];
dgzdata.spk_types    = [dgzdata.spk_types;dgztmp.spk_types];
dgzdata.spk_channels = [dgzdata.spk_channels;dgztmp.spk_channels];
dgzdata.spk_times    = [dgzdata.spk_times;dgztmp.spk_times];
if isfield(dgztmp,'obs_times')
  dgzdata.obs_times    = [dgzdata.obs_times;dgztmp.obs_times];
end

%for i=1:cur_obsnum
%  dgzdata.spk{old_obsnum+i}   = dgztmp.spk{i};
%  dgzdata.spkch{old_obsnum+i} = dgztmp.spkch{i};
%end
