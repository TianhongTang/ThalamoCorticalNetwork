function header = iplx_headeredit(filename)

% Header definition (Laid out for easy reading)
header = struct( ...
    'Filename', filename, ...      % Name of the file from which the header
    ...                            % was read
    'MagicNumber', [], ...         % Always 0x58454c50
    'Version', [], ...             % File format version
    'Comment', '', ...             % User comment, max 128 characters
    'ADFrequency', [], ...         % Timestamp frequency, Hz
    'NumDSPChannels', [], ...      % Number of DSP (spike) channel headers
    'NumEventChannels', [], ...    % Number of event channel headers
    'NumSlowChannels', [], ...     % Number of A/D channel headers
    'NumPointsWave', [], ...       % Number of data points per waveform
    'NumPointsPreThr', [], ...     % Number of data points before threshold
    ...                            % crossing
    'Date', [], ...                % MATLAB "datevec" style, YMDHMS
    'FastRead', [], ...            % Unknown
    'WaveformFreq', [], ...        % Waveform sampling frequency
    'LastTimestamp', [], ...       % Duration of experiment, in ticks
    'Trodalness', [], ...          % 1 for single, 2 for stereotrode, 4 for
    ...                            % tetrode
    'DataTrodalness', [], ...      % Trodalness of the data representation
    'BitsPerSpikeSample', [], ...  % ADC resolution for spike waveforms
    'BitsPerSlowSample', [], ...   % ADC resolution for slow channel data
    'SpikeMaxMagnitudeMV', [], ... % Zero-to-peak of spike waveform (in mV
    ...                            % as ADC values)
    'SlowMaxMagnitudeMV', [], ... % Zero-to-peak of slow data (in mV as ADC
    ...                           % values)
    'SpikePreAmpGain', [], ...
    ... % Channel-by-unit matrices containing the number of spike
    ... % timestamps and waveforms for each channel/unit:
    'TSCounts', [], ... % Number of timestamps
    'WFCounts', [], ... % Number of waveforms
    ... % Starting at event channel 300, EVCounts also includes the number
    ... % of total samples on the continuous channels:
    'EVCounts', [], ... % Number of events on kth event channel
    ... Struct vectors of header information:
    'SpikeHeaders', struct, ... % Spike channel headers
    'EventHeaders', struct, ... % Event channel headers
    'SlowHeaders', struct ...  % Slow channel headers
    );

header.SpikeHeaders = struct( ...
    'Name', '', ...    % DSP Channel name
    'SIGName', '', ... % Corresponding SIG channel name (?)
    'Channel', [], ... % DSP channel number, 1-based
    'WFRate', [], ...  % Limit of waveforms/second divided by 10, when MAP
    ...                % is doing waveform rate limiting
    'SIG', [], ...     % Associated SIG channel number, 1-based
    'Ref', [], ...     % SIG channel used as reference signal, 1-based
    'Gain', [], ...    % Actual gain divided by SpikePreAmpGain (pre
    ...                % version 105, actual gain divided by 1000)
    'Filter', logical([]), ...
    'Threshold', [], ... % Spike detection threshold in AD values
    'Method', [], ...    % Method used to sort units; 1=boxes, 2=templates
    'NUnits', [], ...    % Number of sorted units
    'Template', [], ...  % Template used for template sorting, in AD values
    'Fit', [], ...       % Template fit
    'SortWidth', [], ... % How many points to use in template sorting
    'Boxes', [], ...     % Boxes used in boxes sorting
    'SortBeg', [], ...   % Beginning of sorting window to use in template
    ...                  % sorting
    'Comment', '' ...
    );
header.EventHeaders = struct( ...
    'Name', '', ...    % Event name
    'Channel', [], ... % Event number, 1-based
    'Comment', '' ...
    );
header.SlowHeaders = struct( ...
    'Name', '', ...         % Channel name
    'Channel', [], ...      % Channel number, 0-based
    'ADFreq', [], ...       % Digitization frequency
    'Gain', [], ...         % Gain at the ADC card
    'Enabled', logical([]), ...
    'PreAmpGain', [], ...   % Gain at preamp
    'SpikeChannel', [], ... % Associated spike channel number
    'Comment', '' ...
    );

fid = fopen(filename,'r');

in = fopen(fid,2,'*int32');
header.MagicNumber = typecast(in(1),'uint32');
header.Version = double(in(2));
in = fread(fid,128,'*char');
header.Comment = CStyleStr(in');
in = fread(fid,14,'int32');
header.ADFrequency = in(1);
header.NumDSPChannels = in(2);
header.NumEventChannels = in(3);
header.NumSlowChannels = in(4);
header.NumPointsWave = in(5);
header.NumPointsPreThr = in(6);
header.Date = in(7:12)';
header.FastRead = in(13);
header.WaveformFreq = in(14);
in = fread(fid,1,'double');
header.LastTimestamp = in;
if header.Version >= 103
    in = fread(fid,4,'char');
    header.Trodalness = in(1);
    header.DataTrodalness = in(2);
    header.BitsPerSpikeSample = in(3);
    header.BitsPerSlowSample = in(4);
    in = fread(fid,3,'uint16');
    header.SpikeMaxMagnitudeMV = in(1);
    header.SlowMaxMagnitudeMV = in(2);
    if header.Version >= 105
        header.SpikePreAmpGain = in(3);
    end
end
fseek(fid,256,'bof');
header.TSCounts = fread(fid,[5 130],'int32');
header.WFCounts = fread(fid,[5 130],'int32');
header.EVCounts = fread(fid,512,'int32');

for k=1:header.NumDSPChannels
    in = fread(fid,64,'*char');
    header.SpikeHeaders(k).Name = CStyleStr(in(1:32)');
    header.SpikeHeaders(k).SIGName = CStyleStr(in(33:64)');
    in = fread(fid,9,'int32');
    header.SpikeHeaders(k).Channel = in(1);
    header.SpikeHeaders(k).WFRate = in(2);
    header.SpikeHeaders(k).SIG = in(3);
    header.SpikeHeaders(k).Ref = in(4);
    header.SpikeHeaders(k).Gain = in(5);
    header.SpikeHeaders(k).Filter = logical(in(6));
    header.SpikeHeaders(k).Threshold = in(7);
    header.SpikeHeaders(k).Method = in(8);
    header.SpikeHeaders(k).NUnits = in(9);
    in = fread(fid,[5,64],'int16');
    header.SpikeHeaders(k).Template = in;
    in = fread(fid,6,'int32');
    header.SpikeHeaders(k).Fit = in(1:5);
    header.SpikeHeaders(k).SortWidth = in(6);
    in = fread(fid,prod([5 2 4]),'int16');
    header.SpikeHeaders(k).Boxes = reshape(in,[5,2,4]);
    in = fread(fid,1,'int32');
    header.SpikeHeaders(k).SortBeg = in;
    in = fread(fid,128,'*char');
    header.SpikeHeaders(k).Comment = CStyleStr(in');
    fseek(fid,44,'cof');
end

for k=1:header.NumEventChannels
    in = fread(fid,32,'*char');
    header.EventHeaders(k).Name = CStyleStr(in');
    in = fread(fid,1,'int32');
    header.EventHeaders(k).Channel = in;
    in = fread(fid,128,'*char');
    header.EventHeaders(k).Comment = CStyleStr(in');
    fseek(fid,132,'cof');
end

for k=1:header.NumSlowChannels
    in = fread(fid,32,'*char');
    header.SlowHeaders(k).Name = CStyleStr(in');
    if header.Version >= 104
        numToRead = 6;
        extraPad = 0;
    else
        numToRead = 5;
        extraPad = 4;
    end
    in = fread(fid,numToRead,'int32');
    header.SlowHeaders(k).Channel = in(1);
    header.SlowHeaders(k).ADFreq = in(2);
    header.SlowHeaders(k).Gain = in(3);
    header.SlowHeaders(k).Enabled = logical(in(4));
    header.SlowHeaders(k).PreAmpGain = in(5);
    if header.Version >=104
        header.SlowHeaders(k).SpikeChannel = in(6);
    end
    in = fread(fid,128,'*char');
    header.SlowHeaders(k).Comment = CStyleStr(in');
    fseek(fid,112+extraPad,'cof');
end

fclose(fid);

function str = CStyleStr(str)
endIndex = find(str==char(0),1);
if ~isempty(endIndex)
    if endIndex==1
        str='';
    else
        str = str(1:(endIndex-1));
    end
end