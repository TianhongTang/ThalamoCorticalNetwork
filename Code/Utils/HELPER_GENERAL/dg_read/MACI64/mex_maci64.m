%
% mex_glx64.m : batch file to make all mex DLLs
%
% DAL 16-Sept-07

mex -I. -I/usr/include -I/Applications/MATLAB_R2010b.app/extern/include -L/usr/lib dg_read.c dynio.c df.c dfutils.c flip.c /usr/lib/libz.dylib /Applications/MATLAB_R2010b.app/bin/maci64/libmx.dylib
