function tdat = PVrdFid2(dir, filenum, optin);
% %  dat = PVrdFid(dir, filenum, optin)
% %
% %  Bruker ParaVision:
% %  read FID file
% %
% %  default options:
% %     opt(  'NS',[0 0],       % slice range
% %           'NR',[0 0],       % range in time series
% %           'XY',[0 0],       % single xy point
% %           'VERBOSE','1')
% %
% %  return: series of FID data (in PRECISION type: NOT DOUBLE !!! -> faster)
% %          global acqp/reco parameters
% %
% %  Tested for: onepulse, EPI
% %
% %  Sep 2001 -  Josef Pfeuffer
% %
FCTNAME = 'PVrdFid';

global STDPATH
global acqp

FIDFILE = 'fid';        
             % !! watch out in case of 'ser': 1024 byte minimum blocksize
PLOTRATE = 10;          % data to plot
f_complex = 1;          % reformat re/im data to complex type 

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.NS = [0 0];
dopt.NR = [0 0]; 
dopt.XY = [0 0];
dopt.FUDGE = [0 0 0 0];   % fudged dims for unknown BRUKER data format features
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

nsrange = dopt.NS;
nrrange = dopt.NR; 
xy = dopt.XY;
fudgeDim = [0 0 0 0];
fudgeDim(1:length(dopt.FUDGE)) = dopt.FUDGE;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

if (filenum > 0)    % stdpath handling
   file = sprintf('%s%s\\%s\\%s', STDPATH.pv,num2str(dir),num2str(filenum),FIDFILE);
else              % give full path as arg1
   file = sprintf('%s\\%s', num2str(filenum), FIDFILE);
end
      % read Parameter file(s)
acqp = PVrdPar(dir, filenum);

nx = acqp.ACQ_size(1);
if (acqp.ACQ_dim == 2)
   ny = acqp.ACQ_size(2);
else
   ny = 1;
end
nslices = acqp.NSLICES;
nr = acqp.NR;
%nr = acqp.NI;
if strcmp(acqp.ACQ_word_size, '_16_BIT')
   NBYTE = 2;              % 16 bit 
   PRECISION = 'int16';
elseif strcmp(acqp.ACQ_word_size, '_32_BIT')
   NBYTE = 4;              % 32 bit 
   PRECISION = 'int32';
else
   acqp.ACQ_word_size
   error(sprintf('%s: unknown type <%s>', FCTNAME, acqp.ACQ_word_size))
end
if strcmp(acqp.BYTORDA, 'big')
   BYTEORDER = 'ieee-be';     % 'ieee-be' big endian for PC
elseif strcmp(acqp.BYTORDA, 'little')
   BYTEORDER = 'ieee-le';     % 'ieee-le' little endian for Linux ??
else
   acqp.BYTORDA
   error(sprintf('%s: unknown type <%s>', FCTNAME, acqp.BYTORDA))
end
      % consistency check
fsize = checkfilesize(file);
if fudgeDim(1) > 0
    nx = fudgeDim(1);
end
if fudgeDim(2) > 0
    ny = fudgeDim(2);
end
if fudgeDim(3) > 0
    nslices = fudgeDim(3);
end
if fudgeDim(4) > 0
    nr = fudgeDim(4);
end
if ( fsize/NBYTE == nx*ny*nslices*nr)
   fprintf('reading <%s>\n', file);
   fprintf('%s %s\n', acqp.PULPROG, acqp. GRDPROG);
   fprintf('dim <%d><%d><%d><%d>  filesize %dk\n', ...
            nx, ny, nslices, nr, fsize/1024);
else
   error(sprintf('%s: nint <%d> ~= dimensions <%d><%d><%d><%d>\n', ...
            FCTNAME, fsize/NBYTE, nx, ny, nslices, nr));
end   

fid = fopen(file, 'r', BYTEORDER);
         %%%%%%%%%%%  read all data  %%%%%%%%%%%%%%
if (nsrange(1) == 0 & nrrange(1) == 0 & xy(1) == 0) 
   flen = nx*ny*nslices*nr;
   [tdat, fcount] = fread(fid, flen, [PRECISION '=>' PRECISION]); 
   if (fcount < flen)
      error(sprintf('%s: only %d / %d read', FCTNAME, fcount, flen))
   end
else     %%%%%%%%%%%  read slice / nr range  %%%%%%%%%%%%%%

   if (nsrange(1) ~= 0)
       ns1 = max( [nsrange(1) 1] );
       ns2 = min( [nsrange(2) nslices] );
       if (ns2 < ns1) ns2 = nslices; end
	else
       ns1 = 1;
       ns2 = nslices;
	end
	if (nrrange(1) ~= 0)
       nr1 = max( [nrrange(1) 1] );
       nr2 = min( [nrrange(2) nr] );
       if (nr2 < nr1) nr2 = nr; end
    else
       nr1 = 1;
       nr2 = nr;
	end
	nslices_12 = ns2 - ns1 + 1;
	nr_12 = nr2 - nr1 + 1;
	
   if (xy(1) == 0)
      if (nslices_12 == nslices)
			fseek(fid, NBYTE*nx*ny*nslices*(nr1-1), 'bof');
          flen = nx*ny*nslices*nr_12;
          [tdat, fcount] = fread(fid, flen, [PRECISION '=>' PRECISION]); 
          if (fcount < flen)
             error(sprintf('%s: only %d / %d read', FCTNAME, fcount, flen))
          end
      else
         error('NOT yet Debugged');
                % --- Move to beginning of first time point nr1. 
     		fseek(fid, NBYTE*nx*ny*nslices*(nr1-1), 'bof');
                % --- Move to beginning of first selected slice ns1
         fseek(fid, NBYTE*nx*ny*(ns1-1), 'cof');
	
         flen = nx*ny*nslices_12;
		   for inr=1:nr_12,
             [tdat0, fcount] = fread(fid, flen, [PRECISION '=>' PRECISION]); 
             if (fcount < flen)
                error(sprintf('%s - fread: only %d / %d read', FCTNAME, fcount, flen))
             end
		      if (inr <= 1)
                tdat = tdat0;
             else
                tdat = [tdat tdat0];
                size tdat
             end
                % --- Move to beginning of next nr block:
			   fseek(fid, NBYTE*nx*ny*(nslices-nslices_12), 'cof');
	      end
      end
   else     % read only single FID point (x,y)
      xy1 = min( [nx xy(1)] );
      xy2 = min( [ny xy(2)] );
      if (f_complex)
         npts2read = 2;
      else
         npts2read = 1;
      end

      for inr=nr1:nr2
	    for ins=ns1:ns2
               % --- Move to beginning of time point 
 		   fseek(fid,  NBYTE*nx*ny*nslices*(inr-1) + ...
 		               NBYTE*nx*ny*(ins-1) + ...
 		               NBYTE*nx*(xy2-1) + ...
                     NBYTE*(xy1-1)*npts2read , 'bof');
         [tdat0, fcount] = fread(fid, npts2read, [PRECISION '=>' PRECISION]); 
         if (fcount < npts2read)
            error(sprintf('%s - fread: only %d / %d read', FCTNAME, fcount, flen))
         end
	      if (inr <= 1)
            tdat = zeros(npts2read, nslices_12, nr_12);
         end
         tdat(:,ins-ns1+1,inr-nr1+1) = tdat0;
       end
      end
      nx = npts2read;
      ny = 1;
   end  
   nslices = nslices_12;
   nr = nr_12;
end 
fclose(fid);

if (f_complex)
   fprintf('resorting to COMPLEX ...\n');
   tdat = reshape(tdat,2,nx/2*ny*nslices*nr);
   tdat2 = tdat(2,:);
   tdat = tdat(1,:);
   tdat = complex(tdat, tdat2);
   fprintf('reshaping tdat ...\n');
   tdat = reshape(tdat, nx/2, ny, nslices, nr);
   nx = nx/2;
else
   fprintf('reshaping tdat ...\n');
   tdat = reshape(tdat, nx, ny, nslices, nr);
end
fprintf('dim <%d><%d><%d><%d> \n', nx, ny, nslices, nr);

if f_verbose
figure(3);
colormap(gray);
subplot(1,1,1);
if (ny > 1)
   plotdat = double(tdat(:,:,1,1));
   imagesc(abs(plotdat));
   title('k-space ABS()')
elseif (nx > 1)
   plotdat = squeeze(double(tdat(:,:,1,1)));
   plot(abs(plotdat));   
   hold on
   plot(real(plotdat));   
   hold off
   title('FID abs/real')
else
   plotdat = squeeze(double(tdat(:,:,1,:)));
   plot(abs(plotdat));   
   title('timecourse')
end
end


