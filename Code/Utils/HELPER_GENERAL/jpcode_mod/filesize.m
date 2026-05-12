function fsize = filesize(file);
%
%  determines filesize [byte]
%
%  Sep 2001 -  JP
%

   % --- arg handling
narg = nargin;
error(nargchk(1,1,narg));

fid = fopen(file,'r');
if (fid <= 0)
   error(sprintf('filesize: can not open <%s>', file));
end

status = fseek(fid,0,'eof');
fsize = ftell(fid);     % total bytes