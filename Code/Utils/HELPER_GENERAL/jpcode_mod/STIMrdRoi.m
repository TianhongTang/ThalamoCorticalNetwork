function roi = STIMrdRoi(file)
%
%  STIMULATE:
%%   reads ROI file,  returns coordinates
%%   roi = readroi(file)
%%   
%%   JP Apr 2000

narg = nargin;
error(nargchk(1,1,narg));
		   
LONGLINESTRING = 'period';

fprintf('reading <%s> ...\n', file );
[fid, errmsg] = fopen(file,'r');
if (fid <= 2)
   error(errmsg);
end

boxidx = 0;
while ( ~feof(fid) )
   boxidx = boxidx + 1;
   linearr(boxidx,:) = fgetl(fid);
end

fclose(fid);

if (boxidx > 1)
   linearr
   warning('!! multiple Boxes NOT supported');
elseif (boxidx < 1)
   error('!! no Boxes found');
end

boxstr = linearr(1,:);
fprintf('%s\n', boxstr);
[roi, num] = strsplit(boxstr);
if (num ~= 6)       % e.g. 'Box 3 21 70 94 165'
   error('! NOT a valid Box format !');
end
roi = str2num(char( roi(2:6,:) ));

%--------------------------------------------------------------

