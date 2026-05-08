% IPLX  A small toolbox for the reading of Plexon .plx files.
% This toolbox is designed as a replacement for the m-files and MEX file
% provided by Plexon to read in .plx files. Its primary advantage is that,
% unlike the toolbox provided by Plexon, it is capable of running on
% non-Windows platforms. It is also considerably faster, and provides a
% more "MATLAB-like" interface. The amount of code encapsulated into a MEX
% file has been minimized to solely those components that pure MATLAB is
% too slow at running. To port IPLX to a new platform (in particular if it
% generates errors relating to "iplx_datascan"), make sure MEX is
% configured on your system and then type "mex iplx_datascan.c". This
% should be all you need to do.
% 
% The IPLX toolbox is intended to be more or less compatible with the
% standard Plexon functions. The following functions are provided and
% should have identical interfaces to the Plexon functions of the same
% names (if not, that is a bug).
%    plx_ad
%    plx_ad_v
%    plx_ad_chanmap
%    plx_ad_gap_info
%    plx_ad_span
%    plx_ad_span_v
%    plx_adchan_freqs
%    plx_adchan_gains
%    plx_adchan_names
%    plx_adchan_samplecounts
%    plx_chan_filters
%    plx_chan_gains
%    plx_chan_names
%    plx_chan_thresholds
%    plx_chanmap
% ** plx_close **
%    plx_event_names
%    plx_event_ts
%    plx_info
%    plx_ts
%    plx_waves
%    plx_waves_v
% **Note that it is very important to call PLX_CLOSE after you are done
% reading data from a file. The IPLX toolbox loads the entire contents of
% the Plexon file into persistent memory to mimic the behavior of the
% original plx_* functions, and that memory is not freed until PLX_CLOSE is
% called. The IPLX toolbox only keeps one file in persistent memory at a
% time, so it is much faster to read everything from a file at once rather
% than switch back and forth between files.
% 
% The following functions are not currently supported:
%    plx_information
%    plx_vt_interpret
% If you need an unsupported function, please let me (Geoff) know, and I
% will try to add it.
% 
% The above functions are included for backward compatibility with the old
% Plexon functions. However, when writing new code, it may be preferable to
% use the IPLX specific functions. In particular,
%    iplx_fullfile
% can be used to read the entire file into a single struct for convenient
% access (see IPLX_FULLFILE for more information).
% A few other functions are provided to work on the struct returned by
% IPLX_FULLFILE:
%    iplx_tidychans        Remove channels with no data
%    iplx_multislowchan    Extract all slow channel data
%    iplx_getunit          Extract timestamps & waveforms for a single unit
% There are also low-level functions that typically do not need to be used
% directly if IPLX_FULLFILE is used:
%    iplx_header         Read all file header information into a struct
%    iplx_dataheaders    Read all data block headers into a struct
%    iplx_getspikes
%    iplx_getevents
%    iplx_getslow
%    iplx_datascan       MEX file to extract data blocks
%    iplx_mimichelper    Helper function for the plx_* mimic functions
