% IPLX_DATASCAN  Scan data from a Plexon file to get block information
% Usage:
%    blockLocs = iplx_datascan(datastream)
%    [blockLocs, allData] = iplx_datascan(datastream)
%
% datastream is a vector of int16 values read from a Plexon .plx file. The
% datastream vector should be read starting from the very end of the header
% information (file header + channel headers) to the end of the file.
%
% blockLocs is a vector of indices to the datastream indicating the start
% of each data block header. MATLAB code can then be used to interpret the
% data headers at these locations (see IPLX_DATAHEADERS).
%
% allData is a cell array, with each element corresponding to the data
% contained in one data block.
%
% IPLX_DATASCAN is a MEX, because the nature of PLX files is such that it
% requires fully serial reading that necessitates a for loop; it cannot be
% vectorized. For very large files, attempting to do this in an M-file is
% extraordinarily slow, so it is necessary to use a MEX to parse the data
% blocks. However, once this has been done, MATLAB is perfectly speedy at
% interpreting the header information, etc., so IPLX_DATASCAN does not
% attempt to do this.