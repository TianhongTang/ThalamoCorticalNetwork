function epi_dork( fpath, grfx_on, partialdork_on )
%
% Correction of Dynamic Off-Resonance in K-space (DORK)
%      J. Pfeuffer, P.-F. van de Moortele, X. Hu, G. H. Glover, Magn.Res.Med 2001
%
% Based on rawdork.m by Gary Glover
%  rev xxx	3/30/01	ghg	fix the ref scan, make navrf = center volume
%  rev 1	6/01/01	jp      flag partial dork; TE/TN taken from correct position
%                               option to switch partial DORK; write phiI + phiN
%


if (nargin > 1 )
   grfxBool = 1
else
   grfxBool = 0
end

if (nargin > 2 )
   flag_partialdork = 1      % if != 0:  forces to do only PARTIAL dork correction
   grfxBool = 0
else
   flag_partialdork = 0
end

format compact 
help epi_dork

   % Check for the acquisition timing
atPath = [fpath '/timing_old.float64'];   
if (2 ~= exist( atPath )) error( ['Not found: ', atPath]); end

ppPath = [fpath '/procpar.orig'];
if (2 ~= exist( ppPath )) error( ['Not found: ', ppPath]); end

fidPath = [fpath '/fid.orig'];
if (2 ~= exist( fidPath )) error( ['Not found: ', fidPath]); end

forfixPath = [fpath '/fid.forfix'];
cmd = sprintf( '/bin/cp -f %s %s', fidPath, forfixPath );
if (0 ~= unix( cmd )) error( ['Unix cmd failed:', cmd]); end

cmd = sprintf( '/bin/chmod +w %s', forfixPath );
if (0 ~= unix( cmd )) error( ['Unix cmd failed:', cmd]); end

fp = fopen( forfixPath, 'r+', 'ieee-be');
if (-1 == fp) error( ['Open failed : ', forfixPath]); end

  mh   = getMainHdr(fp)
  np   = getPPV( 'np', ppPath)
  pss  = getPPV( 'pss', ppPath)
  nv   = getPPV( 'nv', ppPath)
  nseg = getPPV( 'nseg', ppPath)
  navecho = getPPV( 'navecho', ppPath )
  te_ppv = getPPV( 'te', ppPath)
  ref_scan_num = getPPV( 'ref_pos', ppPath) + 1;
  fprintf('ref_scan_num = %d\n', ref_scan_num);

  dx = np / 2;
  num_slices = size(pss, 2);

  mainHdrSize = 32;
  blkHdrSize = 28;
  offArr = zeros(1, mh.ntraces * mh.nblocks);
  fileOff = mainHdrSize;
  idx = 1;
  for blkIdx=1:mh.nblocks
     fileOff = fileOff + blkHdrSize;
     for trcIdx=1:mh.ntraces
        offArr(idx) = fileOff;
          idx = idx + 1;
        fileOff = fileOff + mh.tbytes;
     end
  end

  if (mh.ebytes == 2)
     precision = 'int16';
  else
     precision = 'int32';
  end

ntrc = nv/nseg       % Number of traces per segment
nvol = (mh.ntraces * mh.nblocks) / (num_slices * nv);

   % xfer values from JPS to GG variables; laziness here
nx = dx;              % number of complex points
ny = ntrc;            % number of traces per segment
nslc = num_slices;    % number of slices per volume
nint = nseg;          % number of segments per slice
nfr = nvol;           % number of volumes

%wbh = waitbar(0, [num2str(nvol), ' volumes to process']);
% set(wbh, 'Name', 'Correction of Dynamic Off-Resonance in K-space (DORK)' );

fprintf('\ndata: nx ny nseg nslc nvol = %d %d %d %d %d\n\n', ...
	nx, ny, nseg, nslc, nvol);

fprintf('\nReading acqisition timing file: %s\n', atPath);
atfp = fopen( atPath, 'r', 'ieee-be.l64' );
[atArr, atCnt] = fread( atfp, inf, 'float64');
fclose( atfp );

if (atCnt ~= (nx*ny))
   fprintf( 'Read: %d timing points but expected %d = %dx%d\n', ...
	      atCnt, (nx*ny), nx, ny );
   error( ['Wrong number of timing points in: ', atPath]); 
end
atArr = reshape( atArr, nx, ny );    % tells when each acq point was taken

%%%TEplus = 0.000
%%%atArr = atArr + TEplus;

[atmin iatmin] = min(reshape( atArr, 1, nx*ny ));
[atmin2 iatmin2] = min(reshape( atArr(:,2:ny), 1, nx*(ny-1) ));
[atmax iatmax] = max(reshape( atArr, 1, nx*ny ));
atline = atArr(1,3)-atArr(1,2);
fprintf('timing: min min2 max tline = %5.2f %5.2f %5.2f %5.2f ms\n', atmin*1000, atmin2*1000, atmax*1000, atline*1000);

qq = i;
segMtrx = qq( ones(nx,1), ones(ny,1) );

nspv = nseg * nslc;  % Number of segments per volume;
numVolTrc = nspv * ntrc; % Number of traces per volume;
realIdx = 1:2:np;
imagIdx = realIdx + 1;

   % Use a volume near the middle for the reference.
refIdx = round(nvol/2); 

for sIdx=1:nspv
   segOff = ((refIdx-1) * numVolTrc) + ((sIdx-1) * ntrc);

   for tIdx=1:ntrc 
      trcOff = segOff + tIdx;
      if (0 ~= fseek( fp, offArr(trcOff), 'bof'))
        	error('Seek failed'); 
      end
      trc = fread( fp, np, precision );
      segMtrx(:, tIdx) = trc(realIdx) + i*trc(imagIdx);
   end

   [xmax ixmax] = max(abs(segMtrx(:, 2:ntrc)));
   [ymax iepiy] = max(xmax);
   iepix = ixmax(iepiy);
   iepiy = iepiy + 1;   % offset for navecho
   [xmax inecx] = max(abs(segMtrx(:,1)));  % here's nave echo
   fprintf('center of k-space, navecho = (%d %d, %d)\n', iepix, iepiy, inecx)

   iTE = atArr(iepix, iepiy);
   iTN = atArr(iepix, 1);    % don't care about time shift of TN away from middle

   kRef(sIdx).iepix = iepix;
   kRef(sIdx).iepiy = iepiy;
   kRef(sIdx).inecx = inecx;
   kRef(sIdx).epiVal = segMtrx(iepix, iepiy); % get reference epi nav data;
   kRef(sIdx).neVal = segMtrx(inecx, 1);      % get reference epi nav data;
   kRef(sIdx).iTE   = iTE;
   kRef(sIdx).iTN   = iTN;
end

te = median(kRef(sIdx).iTE);    % TAKE same TE for all segments !
tn = median(kRef(sIdx).iTN);
atArr_norm = atArr - te; % Normalize acq time to center of k-space
   
fprintf('TE TN ppTE = (%5.2f %5.2f %5.2f) ms\n', te*1000, tn*1000, te_ppv*1000)

df_arr = zeros(nvol, nspv);
dphi_arr = zeros(nvol, nspv);
dphie_arr = zeros(nvol, nspv);
dphin_arr = zeros(nvol, nspv);

if (1 == grfxBool)
   pfh = figure( 'Name', 'Correction of Dynamic Off-Resonance in K-space (DORK)');
	subplot(2,1,1)
   xlabel('time, image')
   ylabel('delta freq, Hz')
   title('Frequency change');
	subplot(2,1,2)
   xlabel('time, image')
   ylabel('delta phase, Deg')
   title('Phase change');
	drawnow;
end

fprintf( '\n%d Volumes to process: ', nvol );
for vIdx=1:nvol
   fprintf( '.' );
%	waitbar( vIdx / nvol );
   for sIdx=1:nspv
      segOff = ((vIdx - 1) * numVolTrc) + ((sIdx-1) * ntrc);
   		%% Read the data
      for tIdx=1:ntrc 
      	 trcOff = segOff + tIdx;
         if (0 ~= fseek( fp, offArr(trcOff), 'bof'))
            	error('Seek failed'); 
         end;
         trc = fread( fp, np, precision );
      	 segMtrx(:, tIdx) = trc(realIdx) + i*trc(imagIdx);
      end

   		% Calculate the correction
      navepr = kRef(sIdx).epiVal;
      iepix = kRef(sIdx).iepix;
      iepiy = kRef(sIdx).iepiy;
      navep = segMtrx(iepix, iepiy);

	     %% FULL DORK calculation
      if isequal( navecho, 'y') & isequal( flag_partialdork, 0)
	  navecr = kRef(sIdx).neVal;
	  inecx = kRef(sIdx).inecx;
	  navec = segMtrx(inecx, 1);
          dphie = angle(navepr/navep);   % epi echo
          dphin = angle(navecr/navec);   % nav echo
          dphie_arr(vIdx,sIdx) = dphie;
          dphin_arr(vIdx,sIdx) = dphin;

          dw = (dphie - dphin)/(te - tn);
          dphi = (te*dphin - tn*dphie)/(te - tn);
      else        %% PARTIAL DORK calculation: dumb first order algo
          dphi = angle(navepr/navep);    % epi echo
          dphie_arr(vIdx,sIdx) = dphi;
	  dw = dphi/te;
      end

   		% Apply the correction
      if(vIdx ~= ref_scan_num)
         for iy =1:ny
           corr = exp(i*(dphi + dw*atArr_norm(:,iy)));
           segMtrx(:,iy) = corr.*segMtrx(:,iy);
         end
      end
      df_arr(vIdx:nvol,sIdx) = dw/(2*pi);
      dphi_arr(vIdx:nvol,sIdx) = dphi*180/pi;
   
   		% Write the data
      for tIdx=1:ntrc 
      	 trcOff = segOff + tIdx;
         if (0 ~= fseek( fp, offArr(trcOff), 'bof'))
            	error('Seek failed'); 
         end;
   	 trc(realIdx) = real(segMtrx(:,tIdx));
   	 trc(imagIdx) = imag(segMtrx(:,tIdx));
   	 if (length(trc) ~= fwrite( fp, trc, precision ))
  	    error( ['Write failed: ', forfixPath]); 
         end
      end

      if ((1 == grfxBool) & (1 == sIdx) & (0 == mod(vIdx,20)) )
	subplot(2,1,1)
		plot(df_arr(10:nvol, sIdx));
        xlabel('time, image')
        ylabel('delta freq, Hz')
        title('Frequency change');
	subplot(2,1,2)
		plot(dphi_arr(10:nvol, sIdx));
        xlabel('time, image')
        ylabel('delta phase, Deg')
        title('Phase change');
	drawnow;
      end

   end %%%% for_sIdx
end %%%%%%%%%% for_vIdx
fprintf( '\n! DORK correction successfully done !\n-------------------------\n' );

if (1 == grfxBool) 
   close( pfh ); 
end;
   % close( wbh );
fclose(fp);

df = df_arr;
dp = dphi_arr;
dpe = dphie_arr;
dpn = dphin_arr;
save([fpath '/df'],  'df')
save([fpath '/dp'],  'dp')
save([fpath '/dpe'], 'dpe')
save([fpath '/dpn'], 'dpn')
  
%-------------------------------------------------------------------------------
function mainHdr = getMainHdr( fp )
% getMainHdr returns the main header of an fid filePointer

mainHdr.nblocks = fread( fp, 1, 'int32');
mainHdr.ntraces = fread( fp, 1, 'int32');
mainHdr.np = fread( fp, 1, 'int32');
mainHdr.ebytes = fread( fp, 1, 'int32');
mainHdr.tbytes = fread( fp, 1, 'int32');
mainHdr.bbytes = fread( fp, 1, 'int32');
mainHdr.transf = fread( fp, 1, 'int16');
mainHdr.status = fread( fp, 1, 'int16');
mainHdr.spare1 = fread( fp, 1, 'int32');

return

%-------------------------------------------------------------------------------
function vals = getPPV( ppName, ppPath )
% Get the values for parameter name in procpar file path
% Usage: vals = getPPV( ppName, ppPath )

% fn = 'I_t.fid/procpar'

fp = fopen( ppPath, 'r');
if (-1 == fp) error( ['Open failed: ', ppPath]); end

done = 0;
vals = [];

while( done == 0 )
   line = fgetl(fp);
   if (line == -1)
      done = 1;
   else
%      if ~isletter(line)
%         error( 'bad format')
%      else
	if (strcmp(line(1),ppName(1)))
         [name, attr] = strtok(line);
         attr = str2num(attr);

         % Read in the values
         line = fgetl(fp);
         [cnt, parm] = strtok(line);
         cnt = str2num(cnt);

         % REAL_VALS
         if (attr(2) == 1) 
            vals = str2num( parm );
         % STRING_VALS
         else 
            vals = dbl_quote_extract( parm );
            while( size(vals,1) ~= cnt )
               line = fgetl(fp);
               vals = char(vals, dbl_quote_extract( line ) );
            end
         end

         if (strcmp(name, ppName))
            break;
			else
				vals=[];
			end

         % Read in the enums
         enum_line = fgetl(fp);
      end
   end
end

fclose( fp );


%-------------------------------------------------------------------------------
function outStr = dbl_quote_extract( inStr )
%dbl_quote_extract - Extract String from between pair of Double Quotes

dqIdx = findstr(inStr, '"');
dqCnt = size(dqIdx,2);

if ((dqCnt == 0) | (1 == mod(dqCnt,2)))
   error( 'Bad string double quote balance')
end 

for idx = 1:2:dqCnt
   off = [dqIdx(idx) + 1, dqIdx(idx+1) - 1];
   if (idx == 1)
      outStr = inStr(off(1):off(2));
   else
      outStr = char(outStr, inStr(off(1):off(2)));
   end
end
%-------------------------------------------------------------------------------


