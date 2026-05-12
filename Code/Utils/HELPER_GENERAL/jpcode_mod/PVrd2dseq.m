function idat = PVrd2dseq(dir, filenum, optin,IMGFILE);
% %  dat = PVrd2dseq(dir, filenum, optin)
% %
% %  Bruker ParaVision:
% %  read 2dseq RECO file
% %
% %  default options:
% %     opt(  'NS',[0 0],       % slice range
% %           'NR',[0 0],       % range in time series
% %           'XY',[0 0],       % single xy point
% %           'RECO',1,         % number of reco dir
% %           'GETINFO',0,      % flag to retrieve 2dseq INFO, also: global acqp/reco
% %           'VERBOSE','1')
% %
% %  return: series of FID data (in PRECISION type: NOT DOUBLE !!! -> faster)
% %          global acqp/reco parameters
% %
% %  Tested for: onepulse, EPI
% %
% %  Sep 2001 -  Josef Pfeuffer
% %
FCTNAME = 'PVrd2dseq';

global STDPATH
global acqp reco

if ~exist('IMGFILE') | isempty(IMGFILE)
  IMGFILE = '2dseq';        
end

             % !! watch out in case of 'ser': 1024 byte minimum blocksize
PLOTRATE = 10;          % data to plot
f_complex = 0;          % reformat re/im data to complex type 
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.NS = [0 0];
dopt.NR = [0 0]; 
dopt.XY = [0 0];
dopt.RECO = 1;
dopt.GETINFO = 0;
dopt.VERBOSE = 0;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+2,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

nsrange = dopt.NS;
nrrange = dopt.NR; 
xy = dopt.XY;
reconum = dopt.RECO;      
f_getinfo = dopt.GETINFO;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

if (filenum > 0)    % stdpath handling
   file = sprintf('%s/%s/%s/pdata/%s/%s', STDPATH.pv,num2str(dir), ...
      num2str(filenum),num2str(reconum),IMGFILE);
else              % give full path as arg1
   file = sprintf('%s/%s', num2str(filenum), IMGFILE);
end
fprintf('Loading raw file %s\n',file);

% read Parameter file(s)
acqp = PVrdPar(dir, filenum, opt('VERBOSE',f_verbose));
reco = PVrdParReco(dir, filenum, opt('RECO',dopt.RECO,'VERBOSE',f_verbose) );


nx = reco.RECO_size(1);
if (length(reco.RECO_size) < 2)  %% 1D case
    ny = 1;
else
    ny = reco.RECO_size(2);
end

if reco.RECO_transposition
  nx = reco.RECO_size(2);
  ny = reco.RECO_size(1);
end

nz = prod(reco.RECO_size)/nx/ny;
if (length(reco.RECO_size) > 2)  %% general routines depend on a 2D structure
    ny2 = ny;                    %% fake here 3D, remember original 'ny2'
    ny = prod(reco.RECO_size)/nx;
end   
if strcmp(reco.RECO_image_type, 'COMPLEX_IMAGE')
   f_complex = 0;
   nx = nx * 2;
elseif strcmp(reco.RECO_image_type, 'MAGNITUDE_IMAGE')
   f_complex = 0;
else
   %error(sprintf('%s: unknown type <%s>', FCTNAME, RECO_image_type))
end

nslices = acqp.NSLICES;
if strncmp(acqp.PULPROG, '<BLIP_epi',9) 
    nechoes = acqp.NECHOES;   
else
    if acqp.ACQ_rare_factor > 0
        nechoes = acqp.NECHOES/acqp.ACQ_rare_factor;   % don't count echoes used for RARE phase encode
    else
        nechoes = acqp.NECHOES;
    end
end
nechoes = max([1 nechoes]);
nr1 = acqp.NR;
nr2 = acqp.NI;
if (nr1 > 1)
    nr = nr1;
elseif strncmp(acqp.PULPROG, '<BLIP_epi',9)
    nr = nr1;
    %% don't remember in which case this is correct ??
%%    nr = nr2;
%%    nslices = 1;
elseif strncmp(acqp.PULPROG, '<mdeft',5)
    nr = nr2;
else
    nr = 1;
end
if strcmp(reco.RECO_wordtype, '_16BIT_SGN_INT')
   NBYTE = 2;              % 16 bit 
   PRECISION = 'int16';
elseif strcmp(reco.RECO_wordtype, '_32BIT_SGN_INT')
   NBYTE = 4;              % 32 bit 
   PRECISION = 'int32';
else
   error(sprintf('%s: unknown type <%s>', FCTNAME, reco.RECO_wordtype))
end
if strcmp(reco.RECO_byte_order, 'bigEndian')
   BYTEORDER = 'ieee-be';     % 'ieee-be' big endian for PC
elseif strcmp(reco.RECO_byte_order, 'littleEndian')
   BYTEORDER = 'ieee-le';     % 'ieee-le' little endian for Linux ??
else
   reco.RECO_byte_order
   error(sprintf('%s: unknown type <%s>', FCTNAME, reco.RECO_byte_order))
end


      % consistency check
fsize = checkfilesize(file);
if ( fsize/NBYTE == nx*ny*nechoes*nslices*nr)
    if (f_verbose)
        fprintf('reading <%s>\n', file);
        fprintf('%s %s\n', acqp.PULPROG, acqp. GRDPROG);
        fprintf('dim <%d><%d><%d><%d>  filesize %dk\n', ...
            nx, ny, nechoes*nslices, nr, fsize/1024);
    end
else
   error(sprintf('%s: nint <%d> ~= dimensions <%d><%d><%d><%d>\n', ...
            FCTNAME, fsize/NBYTE, nx, ny, nechoes*nslices, nr));
end   

if ( f_getinfo )
    info.file = file;
    info.fsize = fsize;
    info.precision = PRECISION;
    info.byteorder = BYTEORDER;
    info.nx = nx;
    info.ny = ny;
    info.nechoes = nechoes;
    if (length(reco.RECO_size) > 2)
        info.ny = ny2;
    end
    if nz <= 1    
        info.nslices = nslices;
    else    %% 3D experiment
        info.nslices = nslices*nz;
    end        
    info.nr = nr;
    idat = info;    % return 2dseq INFO
    return
end


nslices = nechoes*nslices;   % !!treat multi echoes like multi-slice!!

fid = fopen(file, 'r', BYTEORDER);
         %%%%%%%%%%%  read all data  %%%%%%%%%%%%%%
if (nsrange(1) == 0 & nrrange(1) == 0 & xy(1) == 0) 
   flen = nx*ny*nslices*nr;
   [idat, fcount] = fread(fid, flen, [PRECISION '=>' PRECISION]); 
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
         [idat, fcount] = fread(fid, flen, [PRECISION '=>' PRECISION]); 
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
             [idat0, fcount] = fread(fid, flen, [PRECISION '=>' PRECISION]); 
             if (fcount < flen)
                error(sprintf('fread: only %d / %d read', fcount, flen))
             end
		      if (inr <= 1)
                idat = idat0;
             else
                idat = [idat idat0];
                size idat
             end
                % --- Move to beginning of next nr block:
			   fseek(fid, NBYTE*nx*ny*(nslices-nslices_12), 'cof');
	      end
      end
   else     % read only single FID point (x,y)
      xy1 = min( [nx xy(1)] );
      xy2 = min( [ny xy(2)] );
      if (f_complex) npts2read = 2; else npts2read = 1; end
      
	   for inr=nr1:nr2
	    for ins=ns1:ns2
               % --- Move to beginning of time point 
 		   fseek(fid,  NBYTE*nx*ny*nslices*(inr-1) + ...
 		               NBYTE*nx*ny*(ins-1) + ...
 		               NBYTE*nx*(xy2-1) + ...
                     NBYTE*(xy1-1)*npts2read , 'bof');
         [idat0, fcount] = fread(fid, npts2read, [PRECISION '=>' PRECISION]); 
         if (fcount < npts2read)
            error(sprintf('fread: only %d / %d read', fcount, flen))
         end
	      if (inr == nr1)
            idat = zeros(npts2read, nslices_12,nr_12);
         end
         idat(:,ins-ns1+1,inr-nr1+1) = idat0;
	   end
     end
      nx = npts2read;
      ny = 1;
   end  
   nslices = nslices_12;
   nr = nr_12;
end 
fclose(fid);

if (length(reco.RECO_size) > 2)
    ny = ny2;
end
if nz > 1    %% 3D experiment
    nslices = nslices*nz;
end

if (f_complex)
   if f_verbose fprintf('resorting to COMPLEX ...\n'); end
   idat = reshape(idat,2,nx/2*ny*nslices*nr);
   idat2 = idat(2,:);
   idat = idat(1,:);
   idat = complex(idat, idat2);
   if f_verbose fprintf('reshaping idat ...\n'); end
   idat = reshape(idat, nx/2, ny, nslices, nr);
   nx = nx/2;
elseif strcmp(reco.RECO_image_type, 'COMPLEX_IMAGE')
   if f_verbose fprintf('resorting to COMPLEX_IMAGE ...\n'); end
   idat = reshape(idat,nx/2*ny*nslices*nr,2);
   idat2 = idat(:,2);
   idat = idat(:,1);
   idat = complex(idat, idat2);
   if f_verbose fprintf('reshaping idat ...\n'); end
   idat = reshape(idat, nx/2, ny, nslices, nr);
   nx = nx/2;
else
   if f_verbose fprintf('reshaping idat ...\n'); end
   idat = reshape(idat, nx, ny, nslices, nr);
end
fprintf('dim <%d><%d><%d><%d>\n', nx, ny, nslices, nr);

if (f_verbose)
	figure(3);
	colormap(gray);
	subplot(1,1,1);
	if (ny > 1)
       ploidat = double(idat(:,:,1,1));
       imagesc(abs(ploidat));
       title('2dseq ABS()')
	elseif (nx > 1)
       ploidat = squeeze(double(idat(:,:,1,1)));
       plot(abs(ploidat));   
       hold on
       plot(real(ploidat));   
       hold off
       title('FID abs/real')
	else
       ploidat = squeeze(double(idat(:,:,1,:)));
       plot(abs(ploidat));   
       title('timecourse')
	end
end


