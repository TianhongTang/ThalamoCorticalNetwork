function fullplx = iplx_tidychans(fullplx)

hasspikes = ~cellfun(@isempty,{fullplx.Spikes.Timestamps});
hasevents = ~cellfun(@isempty,{fullplx.Events.Timestamps});
hasslow   = ~cellfun(@isempty,{fullplx.Slow.Timestamps});
fullplx.Spikes = fullplx.Spikes(hasspikes);
fullplx.Events = fullplx.Events(hasevents);
fullplx.Slow   = fullplx.Slow(hasslow);