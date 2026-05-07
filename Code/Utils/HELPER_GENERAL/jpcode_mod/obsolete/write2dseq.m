function flag=Write2dseq(inarray,fname,precision,byteorder);
%
% flag=Write2dseq(inarray,fname,precision,byteorder)
% inarray			: n-dimensional Matlab array to write
% fname		: name of output file
% precision		: data type specifier
% byteorder		: (s) swap, (n) non-swap is required
%
% mWrite2dseq creates an unformatted binary file in the current
% current working directory. Used to create output similar to
% Bruker 2dseq files.

if (byteorder == 's'),
	fid=fopen(fname,'wb','ieee-be');
else
	fid=fopen(fname,'wb');
end

% get sturcture of array:
dims = length(size(inarray));

if (dims < 3),
	error('mWrite2dseq: Expecting a 3D or 4D input array');
	return;
end

% number of elements in array

ne = prod(size(inarray));

% transpose array and reshape
%inarray = farray';
inarray = squeeze(reshape(inarray,1,ne));

count = fwrite(fid,inarray,precision);

if ((ne - count) == 0)
   flag = 1;
   fprintf('Wrote binary file: %s\n',fname);
else
   flag = 0;
   fprintf('mWrite2dseq could not write: %s\n',fname);
end


fclose(fid);
