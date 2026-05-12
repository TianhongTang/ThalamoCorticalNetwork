function eventdata = iplx_getevents(hdr,data)

% eventdata fields:
%    Name
%    Channel
%    Comment
%    Events
%        TimeStamp
%        Strobed

eventdata = hdr.EventHeaders;

data = pickFromData(data,data.Type==4);

for k=1:numel(eventdata)
    chdata = pickFromData(data,data.Channel==eventdata(k).Channel);
    if ~isempty(chdata.TimeStamp)
        eventdata(k).Events.TimeStamp = double(chdata.TimeStamp) / ...
            double(hdr.ADFrequency);
        if eventdata(k).Channel==257
            eventdata(k).Events.Strobed   = chdata.Unit;
        else
            eventdata(k).Events.Strobed   = [];
        end
    else
        eventdata(k).Events = struct('TimeStamp',{},'Strobed',{});
    end
end

function newdata = pickFromData(data,pick)
newdata.TimeStamp = data.TimeStamp(pick);
newdata.Channel = data.Channel(pick);
newdata.Unit = data.Unit(pick);