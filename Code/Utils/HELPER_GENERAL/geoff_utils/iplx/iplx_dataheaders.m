function [datahdrs,wordStream] = iplx_dataheaders(header)

datahdrs = struct( ...
    'DataPosition', [], ... % Position of this header in the file
    'Type', uint8([]), ... % 1=spike, 4=event, 5=slow
    'TimeStamp', [], ... % Timestamp in ticks. This is represented in the
    ...                  % PLX file as a 40 bit value; we will convert it
    ...                  % to a double
    'Channel', int16([]), ...
    'Unit', int16([]), ...
    'NumberOfWaveforms', int16([]), ... % Number of waveforms represented in this
    ...                          % data block, usually 0 or 1.
    'NumberOfWordsInWaveform', int16([]), ... % Number of samples per waveform in
    ...                                % this data block.
    'Data',[] ...
    );

fileHeaderBytes = 7504;
spikeHeaderBytes = 1020;
eventHeaderBytes = 296;
slowHeaderBytes = 296;

numSpikeHdrs = header.NumDSPChannels;
numEventHdrs = header.NumEventChannels;
numSlowHdrs = header.NumSlowChannels;
totalHdrSize = fileHeaderBytes + ...
               numSpikeHdrs*spikeHeaderBytes + ...
               numEventHdrs*eventHeaderBytes + ...
               numSlowHdrs*slowHeaderBytes;

fid = fopen(header.Filename,'r');
fseek(fid,totalHdrSize,'bof');

wordStream = fread(fid,'*int16');
fclose(fid);

[dataLocs,data] = iplx_datascan(wordStream);
datahdrs.Data = data;
datahdrs.DataPosition = dataLocs;
datahdrs.Type = uint8(wordStream(dataLocs));
tsWord2 = reshape(wordStream(dataLocs+1),1,[]);
tsWord0 = reshape(wordStream(dataLocs+2),1,[]);
tsWord1 = reshape(wordStream(dataLocs+3),1,[]);
tsWord3 = zeros(size(tsWord0),'int16');
tsMat = [tsWord0;tsWord1;tsWord2;tsWord3];
datahdrs.TimeStamp = double(typecast(tsMat(:),'uint64'));
datahdrs.Channel = wordStream(dataLocs+4);
datahdrs.Unit = wordStream(dataLocs+5);
datahdrs.NumberOfWaveforms = wordStream(dataLocs+6);
datahdrs.NumberOfWordsInWaveform = wordStream(dataLocs+7);
