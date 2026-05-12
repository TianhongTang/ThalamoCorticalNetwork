function dgz = dg_read_local(filename)
% DG_READ_LOCAL  Copy a DGZ file locally in order to open it
% This is a workaround for the strange bug that causes dg_read to produce a
% "file not found" error when reading a DGZ file from Einstein. For
% mysterious reasons, dg_read seems to work if the file is copied locally,
% and dg_read is called from that local directory. DG_READ_LOCAL forces the
% issue by copying the DGZ file to the local temp directory, changing to
% that directory, opening the file, and then restoring MATLAB to its
% original state.

dgreadFile = which('dg_read');
if isempty(dgreadFile)
    error('dg_read not found.');
end
dgreadDir = fileparts(dgreadFile);
origPath = path;
origDir  = pwd;
[dummy, file, ext] = fileparts(filename);
tmpFile = [file ext];
copyfile(filename, tempdir);
cd(tempdir);
addpath(dgreadDir);
try
    dgz = dg_read(tmpFile);
    delete(tmpFile);
    cd(origDir);
    path(origPath);
catch
    delete(tmpFile);
    cd(origDir);
    path(origPath);
end
