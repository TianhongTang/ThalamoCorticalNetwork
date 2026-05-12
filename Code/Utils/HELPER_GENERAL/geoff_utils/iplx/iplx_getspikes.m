function spikedata = iplx_getspikes(hdr,data,channel)

% spikedata fields:
%    Name
%    SIGName
%    Channel
%    WFRate
%    SIG
%    Ref
%    Gain
%    Filter
%    Threshold
%    Method
%    NUnits
%    Template
%    Fit
%    SortWidth
%    Boxes
%    SortBeg
%    Comment
%    Spikes
%        TimeStamp
%        Unit
%        Waveform

if nargin < 3 || isempty(channel)
    spikedata = hdr.SpikeHeaders;
else
    chgood = ismember([hdr.SpikeHeaders.Channel],channel);
    spikedata = hdr.SpikeHeaders(chgood);
end

data = pickFromData(data,data.Type==1);

for k=1:numel(spikedata)
    % Compute the conversion from A/D values to mV:
    switch hdr.Version
        case {100, 101, 102}
            ADConvert = 3000 / (2048 * spikedata(k).Gain * 1000);
        case {103, 104}
            ADConvert = hdr.SpikeMaxMagnitudeMV / ...
                (.5*(2^hdr.BitsPerSpikeSample) * spikedata(k).Gain * 1000);
        otherwise
            ADConvert = hdr.SpikeMaxMagnitudeMV / ...
                ( .5*(2^hdr.BitsPerSpikeSample) * spikedata(k).Gain * ...
                  hdr.SpikePreAmpGain );
    end
    spikedata(k).ADConversion = ADConvert;
    chdata = pickFromData(data,data.Channel==spikedata(k).Channel);
    chdata.Data = cellfun(@(x)(ADConvert*double(x)), chdata.Data, ...
        'UniformOutput',false);
    if ~isempty(chdata.TimeStamp)
        spikedata(k).Spikes.TimeStamp = double(chdata.TimeStamp) / ...
            double(hdr.ADFrequency);
        spikedata(k).Spikes.Unit = chdata.Unit;
        spikedata(k).Spikes.Waveform = [chdata.Data{:}];
    else
        spikedata(k).Spikes = ...
            struct('TimeStamp',{},'Unit',{},'Waveform',{});
    end
end

function newdata = pickFromData(data,pick)
newdata.TimeStamp = data.TimeStamp(pick);
newdata.Channel = data.Channel(pick);
newdata.Unit = data.Unit(pick);
newdata.Data = data.Data(pick);