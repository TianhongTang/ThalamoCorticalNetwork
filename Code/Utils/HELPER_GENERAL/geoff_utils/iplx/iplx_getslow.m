function slowdata = iplx_getslow(hdr,data,channel)

% slowdata fields:
%    Name
%    Channel
%    ADFreq
%    Gain
%    Enabled
%    PreAmpGain
%    SpikeChannel
%    Comment
%    Data
%        TimeStamp
%        Data

if nargin < 3 || isempty(channel)
    slowdata = hdr.SlowHeaders;
else
    chgood = ismember([hdr.SlowHeaders.Channel],channel);
    slowdata = hdr.SlowHeaders(chgood);
end
% slowdata fields:
%   'Name', '', ...         % Channel name
%   'Channel', [], ...      % Channel number, 0-based
%   'ADFreq', [], ...       % Digitization frequency
%   'Gain', [], ...         % Gain at the ADC card
%   'Enabled', logical([]), ...
%   'PreAmpGain', [], ...   % Gain at preamp
%   'SpikeChannel', [], ... % Associated spike channel number
%   'Comment', '' ...
%   'Data', struct

data = pickFromData(data,data.Type==5);
% data fields:
%   'TimeStamp', [], ... % Timestamp in ticks.
%   'Channel', int16([]), ...
%   'NumberOfWaveforms', int16([]), ... % Number of waveforms represented in this
%   ...                          % data block, usually 0 or 1.
%   'NumberOfWordsInWaveform', int16([]), ... % Number of samples per waveform in
%   ...                                % this data block.
%   'Data',[] ...

for k=1:numel(slowdata)
    % Compute the conversion from A/D values to mV:
    switch hdr.Version
        case {100, 101}
            ADConvert=5000 / ( 2048 * slowdata(k).Gain * 1000);
        case 102
            ADConvert=5000 / ( 2048 * slowdata(k).Gain * ...
                slowdata(k).PreAmpGain );
        otherwise
            ADConvert=hdr.SlowMaxMagnitudeMV / ...
               ( .5*(2^hdr.BitsPerSlowSample) * slowdata(k).Gain * ...
                 slowdata(k).PreAmpGain );
    end
    slowdata(k).ADConversion = ADConvert;
    ticksPerSample = hdr.ADFrequency/slowdata(k).ADFreq;
    ch = slowdata(k).Channel;
    chdata = pickFromData(data,data.Channel==ch);
    chdata = sortTime(chdata);
    sampleCounts = chdata.NumberOfWaveforms.*chdata.NumberOfWordsInWaveform;
    tickCounts = sampleCounts(1:end-1)*ticksPerSample;
    deltaTicks = diff(chdata.TimeStamp);
    blockEnds = find(tickCounts~=deltaTicks);
    blockEnds = [blockEnds; numel(chdata.Data)];
    if any(blockEnds)
        for m=1:numel(blockEnds)
            if m == 1
                blockStart = 1;
            else
                blockStart = blockEnds(m-1)+1;
            end
            slowdata(k).Data(m).TimeStamp = ...
                double(chdata.TimeStamp(blockStart)) / ...
                double(hdr.ADFrequency);
            slowdata(k).Data(m).Data = ADConvert * ...
                double(cell2mat(chdata.Data(blockStart:blockEnds(m))));
        end
    else
        slowdata(k).Data = struct('TimeStamp',{},'Data',{});
    end
end

function newdata = pickFromData(data,pick)
newdata.TimeStamp = data.TimeStamp(pick);
newdata.Channel = data.Channel(pick);
newdata.NumberOfWaveforms = data.NumberOfWaveforms(pick);
newdata.NumberOfWordsInWaveform = data.NumberOfWordsInWaveform(pick);
newdata.Data = data.Data(pick);

function data = sortTime(data)
[data.TimeStamp, idx] = sort(data.TimeStamp);
data.Channel = data.Channel(idx);
data.NumberOfWaveforms = data.NumberOfWaveforms(idx);
data.NumberOfWordsInWaveform = data.NumberOfWordsInWaveform(idx);
data.Data = data.Data(idx);