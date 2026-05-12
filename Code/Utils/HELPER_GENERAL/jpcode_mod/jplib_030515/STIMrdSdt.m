function data = STIMrdSdt( fileNameRoot )
%
%  STIMULATE:
% --- writes data in SDT/SPR format  ---
%
% INPUT:
%	data
%	file
%
% OUTPUT:
%	written files
%
% FLAGS:
%	ROTATE: file info
%	VERBOSE: file info
%
%  JP June 2000
%

format compact

BYTEORDER = 'ieee-be';     % 'ieee-be' big endian for PC
%BYTEORDER = 'ieee-le';     % 'ieee-le' little endian for Linux ??

sprPath = [fileNameRoot, '.spr'];
fp = fopen( sprPath ); 
if (fp > 0)
   fprintf('--- reading <%s>\n', sprPath);
else
   error(sprintf('--- can not open <%s>', sprPath));
end
done=0;
while( done == 0 )
   line = fgetl(fp);
   if (line == -1)
      done = 1;
   else
      [attr, val] = strtok(line,':');
      switch attr
         case 'numDim'
            numDimStr = strtok(val,':');
         case 'dim'
            dimStr = strtok(val,':');
         case 'dataType'
            dataTypeStr = strtok(val,':');
      end
   end
end
fclose( fp );

numDimStr
dimStr
dataTypeStr

switch deblank(strjust( dataTypeStr ,'left'))
   case 'REAL'
		precision = 'float32';
      NBYTE = 4;
   case 'WORD'
		precision = 'int16';
      NBYTE = 2;
	otherwise
		error( 'Can not handle dataType')
end

sdtPath = [fileNameRoot, '.sdt'];
dimArr = str2num(dimStr);
      % consistency check
fsize = checkfilesize(sdtPath);
if ( fsize/NBYTE ~= prod(dimArr))
   error(sprintf('%s: nint <%d> ~= dimensions\n', FCTNAME, fsize/NBYTE));
end   
fp = fopen(sdtPath, 'r', BYTEORDER);
%%%dimArr(4) = dimArr(4)-4000;
data = fread( fp, prod(dimArr), precision);
data = reshape( data, dimArr );     % ' rotate image for MATLAB
fclose( fp );
