function fullplx = iplx_fullfile(filename,verbose)
% IPLX_FULLFILE  Read in a full .plx file using the improved PLX library
%
%    fullplx = iplx_fullfile(filename)
%    fullplx = iplx_fullfile(filename,verboseFlag)
% IPLX_FULLFILE reads in the entire .plx file "filename" into the struct
% fullplx. If provided, verboseFlag determines whether or not IPLX_FULLFILE
% should print information about the file as it reads it (default is
% false).
%
% fullplx is a struct with the following fields:
%    Filename
%    Version           Version number (100-105)
%    Date              Time at which the file was recorded
%    FileDuration      Length of the file in seconds
%    Comment           String, max 128 characters
%    TSFreq            Full A/D frequency; number of timestamps/second
%    NumPointsPreThr   Number of pre-threshold samples in spike waveform
%    WaveformFreq      A/D frequency for spike waveforms
%    Trodalness        1=single, 2=stereotrode, 4=tetrode
%    DataTrodalness    
%    Spikes             (# of spike channels)-by-1 struct vector of spike
%                       channel information, with fields:
%        Name            Name of DSP channel
%        Channel         DSP channel #, 1-based
%        SIGName         Name of SIG channel associated with DSP channel
%        SIGChannel      SIG channel # associated with DSP channel, 1-based
%        RefChannel      SIG reference channel, 1-based
%        WFRateLimit     When waveform rate limiting, (limit#/sec)/10
%        ADConversion    Conversion factor from AD values to mV
%        NUnits          Number of sorted units on this channel
%        Filter          Boolean
%        Threshold       Spike detection threshold in mV
%        Sort            Struct of spike sorting information, with fields:
%            Method               1=boxes, 2=templates
%            Template             Templates used, in mV
%            TemplateFit          Eponymous
%            TemplateWidth        # points to use in template sorting
%            TemplateBeginning    Beginning of sorting window
%            Boxes                Boxes used
%        Comment         Spike channel comment
%        Timestamps      Vector of spike timestamps on this channel
%        Units           Corresponding sorted units of the spikes
%        Waveforms       Corresponding waveforms of the spikes,
%                        (# of samples per waveform)-by-(# of spikes)
%    Events            (# of event channels)-by-1 struct vector of event
%                      channel information, with fields:
%        Name          Event name
%        Channel       Event number
%        Comment       Event channel comment
%        Timestamps    Vector of timestamps for this event
%        Strobed       If channel==257, vector of strobed value, else empty
%    Slow              (# of continuous channels)-by-1 struct vector of slow
%                      (continuous) channel information, with fields:
%        Name            Channel name
%        Channel         0-based channel number
%        Enabled         Boolean
%        SampleFreq      Sampling frequency
%        ADConversion    Conversion factor from AD values to mV
%        SpikeChannel    Associated spike channel number
%        Comment         Slow channel comment
%        Timestamps      Vector of contiguous slow data block timestamps
%        Data            Cell vector of slow data, each element is one
%                        contiguous slow data block.
% 
% The large single struct format may not be the most immediately convenient
% for your specific application, but in many cases the data can be
% collected with a few quick commands. Some examples:
%   Get the spike timestamps on all channels as a cell array of timestamp
%   vectors:
%     {fullplx.Spikes.Timestamps}
%   Get the spike waveforms on channel 3, unit 1:
%     isUnit = fullplx.Spikes(3).Units==1;
%     spikes = fullplx.Spikes(3).Waveforms(:,isUnit);

if nargin < 2
    verbose = false;
end
if verbose
    disp('Reading header...');
end
hdr = iplx_header(filename);
if verbose
    disp('Reading all data...');
end
data = iplx_dataheaders(hdr);
if verbose
    disp('Extracting spike information...');
end
spikes = iplx_getspikes(hdr,data);
if verbose
    disp('Extracting event information...');
end
events = iplx_getevents(hdr,data);
if verbose
    disp('Extracting slow waveform data...');
end
slow   = iplx_getslow(hdr,data);

fullplx.Filename        = filename;
fullplx.Version         = hdr.Version;
fullplx.Date            = hdr.Date;
fullplx.FileDuration    = hdr.LastTimestamp / hdr.ADFrequency;
fullplx.Comment         = hdr.Comment;
fullplx.TSFreq          = hdr.ADFrequency;
fullplx.NumPointsPreThr = hdr.NumPointsPreThr;
fullplx.WaveformFreq    = hdr.WaveformFreq;
fullplx.Trodalness      = hdr.Trodalness;
fullplx.DataTrodalness  = hdr.DataTrodalness;

for k=1:numel(spikes);
    ADConv = spikes(k).ADConversion;
    fullplx.Spikes(k).Name         = spikes(k).Name;
    fullplx.Spikes(k).Channel      = spikes(k).Channel;
    fullplx.Spikes(k).SIGName      = spikes(k).SIGName;
    fullplx.Spikes(k).SIGChannel   = spikes(k).SIG;
    fullplx.Spikes(k).RefChannel   = spikes(k).Ref;
    fullplx.Spikes(k).WFRateLimit  = spikes(k).WFRate;
    fullplx.Spikes(k).ADConversion = ADConv;
    fullplx.Spikes(k).NUnits       = spikes(k).NUnits;
    fullplx.Spikes(k).Filter       = spikes(k).Filter;
    fullplx.Spikes(k).Threshold    = spikes(k).Threshold*ADConv;
    fullplx.Spikes(k).Sort.Method            = spikes(k).Method;
    fullplx.Spikes(k).Sort.Template          = spikes(k).Template*ADConv;
    fullplx.Spikes(k).Sort.TemplateFit       = spikes(k).Fit;
    fullplx.Spikes(k).Sort.TemplateWidth     = spikes(k).SortWidth;
    fullplx.Spikes(k).Sort.TemplateBeginning = spikes(k).SortBeg;
    fullplx.Spikes(k).Sort.Boxes             = spikes(k).Boxes;
    fullplx.Spikes(k).Comment      = spikes(k).Comment;
    fullplx.Spikes(k).Timestamps   = spikes(k).Spikes.TimeStamp;
    fullplx.Spikes(k).Units        = spikes(k).Spikes.Unit;
    fullplx.Spikes(k).Waveforms    = spikes(k).Spikes.Waveform;
end

for k=1:numel(events)
    fullplx.Events(k).Name       = events(k).Name;
    fullplx.Events(k).Channel    = events(k).Channel;
    fullplx.Events(k).Comment    = events(k).Comment;
    fullplx.Events(k).Timestamps = events(k).Events.TimeStamp;
    fullplx.Events(k).Strobed    = events(k).Events.Strobed;
end

for k=1:numel(slow)
    fullplx.Slow(k).Name         = slow(k).Name;
    fullplx.Slow(k).Channel      = slow(k).Channel;
    fullplx.Slow(k).Enabled      = slow(k).Enabled;
    fullplx.Slow(k).SampleFreq   = slow(k).ADFreq;
    fullplx.Slow(k).ADConversion = slow(k).ADConversion;
    fullplx.Slow(k).SpikeChannel = slow(k).SpikeChannel;
    fullplx.Slow(k).Comment      = slow(k).Comment;
    fullplx.Slow(k).Timestamps   = [slow(k).Data.TimeStamp];
    fullplx.Slow(k).Data         = {slow(k).Data.Data};
end

if verbose
    fprintf('File: %s (version %d)\n',fullplx.Filename,fullplx.Version);
    fprintf('%s (%g minutes)\n', datestr(fullplx.Date), ...
        fullplx.FileDuration/60);
    fprintf('%s\n', deblank(fullplx.Comment));
    fprintf('Sampling frequency: %d\n',fullplx.TSFreq);
    fprintf('Leading samples in spikes: %d\n',fullplx.NumPointsPreThr);
    fprintf('Spike sampling frequency: %d\n', fullplx.WaveformFreq);
    %disp('Spike channel information:')
    %for k=1:numel(fullplx.Spikes)
    %    
    %    disp(sprintf('    Channel %2d: %d spikes on %d units', ...
    %        fullplx.Spikes(k).Channel, ...
    %        numel(fullplx.Spikes(k).Spikes.TimeStamp), ...
    %        fullplx.Spikes(k).NUnits));
    %end
    
end