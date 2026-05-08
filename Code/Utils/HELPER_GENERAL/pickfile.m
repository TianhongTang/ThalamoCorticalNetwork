function [filename, filedir] = pickfile(msg,fdir,fpatt)
% PURPOSE : To pick up a file
% USAGE : [filename,filedir] = pickfile(msg,fdir,fpatt)
%         'filename' is empty when canceled.
% VERSION : 1.00  18-Sep-2000  YM

if ~exist('msg'), msg = 'Select a file'; end
if ~exist('fdir'), fdir = pwd; end
if ~exist('fpatt'), fpatt = '*.*'; end

filename = [];  filedir = [];
fpatt = sprintf('%s/%s',fdir,fpatt);
[f, d] = uigetfile(strrep(fpatt,'/','\'),msg);
if length(f) == 1 & length(d) == 1, return; end

filename = f;
filedir = strrep(d(1:(length(d)-1)),'\','/');

if nargout == 1,
  filename = sprintf('%s/%s',filedir,filename);
end
