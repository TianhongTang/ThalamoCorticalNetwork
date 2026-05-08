function dgzraw_concat = dgzCat(dgzdata_raw)
%
% Concatenate several data file structures into matlab format
% DAL JAN-2000
%
sz = size(dgzdata_raw);
len = sz(2);

e_pre        = [];
e_times      = [];
e_params     = [];
e_types      = [];
e_subtypes   = [];
ems          = [];
spk_types    = [];
spk_channels = [];
spk_times    = [];

for i=1:len
   e_pre         = [e_pre; dgzdata_raw(i).e_pre];
   e_times       = [e_times; dgzdata_raw(i).e_times];
   e_types       = [e_types; dgzdata_raw(i).e_types];
   e_subtypes    = [e_subtypes; dgzdata_raw(i).e_subtypes];
   e_params      = [e_params; dgzdata_raw(i).e_params];
   ems           = [ems; dgzdata_raw(i).ems];
   spk_types     = [spk_types; dgzdata_raw(i).spk_types];
   spk_channels  = [spk_channels; dgzdata_raw(i).spk_channels];
   spk_times     = [spk_times; dgzdata_raw(i).spk_times];
end

dgzraw_concat.e_pre         = e_pre;
dgzraw_concat.e_times       = e_times;
dgzraw_concat.e_types       = e_types;
dgzraw_concat.e_subtypes    = e_subtypes;
dgzraw_concat.e_params      = e_params;
dgzraw_concat.ems           = ems;
dgzraw_concat.spk_types     = spk_types;
dgzraw_concat.spk_channels  = spk_channels;
dgzraw_concat.spk_times     = spk_times;