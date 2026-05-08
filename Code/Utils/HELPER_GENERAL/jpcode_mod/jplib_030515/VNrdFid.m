function fid = VNrdFid( dn )
%
%  Varian VNMR
%  reads FID data
%

ppPath = [dn '/procpar'];
fidPath = [dn '/fid'];

fp = fopen( fidPath, 'r', 'ieee-be');
fprintf(2, 'Reading: %s\n', fidPath);
mh = VNgetMainHdr( fp );

ni = VNgetPPV( 'ni', ppPath);
np = VNgetPPV( 'np', ppPath);
pss = VNgetPPV( 'pss', ppPath);

dx = np / 2;
dy = ni;
dz = size( pss, 2 );
dt = mh.ntraces * mh.nblocks / (dz * ni);
fprintf( 1, 'Dims: %d x %d x %d x %d\n', dx, dy, dz, dt);

fid = zeros(dx,dy,dz,dt) + i*ones(dx,dy,dz,dt);

[posn, slcOff] = sort( pss );
% fftScale = dx * dy;

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
	precision = 'int16'
else
	precision = 'int32'
end

realIdx = 1:2:np;
imagIdx = realIdx + 1;

for it=1:dt
	fprintf( 2, '\nTime: %d ', it);
	for iz=1:dz
	   fprintf( 2, '.');
		for iy=1:dy
			trcIdx = ((iy-1) * dt * dz) + slcOff(iz) + ((it-1) * dz);
			if (0 ~= fseek( fp, offArr(trcIdx), 'bof'))
				error( 'Seek failed' );
			end
			trc = fread( fp, np, precision );
			fid(:,iy,iz,it) = trc(realIdx) + i*trc(imagIdx);
		end
	end
end
fprintf( 2, '\nDone\n' );

fclose(fp);

%----------------------------------------------------------------------
function blkHdr = VNgetBlkHdr( fp )
%
%  Varian VNMR
%  getBlkHdr -> returns the block header of an fid filePointer

blkHdr.scale  = fread( fp, 1, 'int16');
blkHdr.status  = fread( fp, 1, 'int16');
blkHdr.index  = fread( fp, 1, 'int16');
blkHdr.spare3  = fread( fp, 1, 'int16');
blkHdr.ctcount  = fread( fp, 1, 'int32');
blkHdr.lpval  = fread( fp, 1, 'float32');
blkHdr.rpval  = fread( fp, 1, 'float32');
blkHdr.lvl  = fread( fp, 1, 'float32');
blkHdr.tlt  = fread( fp, 1, 'float32');

return
%----------------------------------------------------------------------
function mainHdr = VNgetMainHdr( fp ) 
%
%  Varian VNMR
%  getMainHdr returns the main header of an fid filePointer

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
%----------------------------------------------------------------------
function vals = VNgetPPV( ppName, ppPath )
%
%  Varian VNMR
%  Get the values for parameter name in procpar file path
%  Usage: vals = getPPV( ppName, ppPath )

% fn = 'I_t.fid/procpar'

fp = fopen( ppPath, 'r');
done = 0;
vals = [];

while( done == 0 )
   line = fgetl(fp);
   if (line == -1)
      done = 1;
   else
      if ~isletter(line)
         error( 'bad format')
      else
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

%----------------------------------------------------------------------
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

%----------------------------------------------------------------------



