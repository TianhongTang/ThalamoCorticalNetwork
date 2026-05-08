function im=Read2dseq(fname,nx,ny,tns,ns1,ns2,nt1,nt2,byteorder);
%
% im=Read2dseq(fname,nx,ny,tns,ns1,ns2,nt1,nt2,byteorder)
% nx,ny		: x- and y- dimensions of image
% tns			: total number of slices in dataset
% ns1,ns2	: slice range
% nt1,nt2	: time point range
% byteorder	: (s) swap, (n) non-swap is required

fname

if (byteorder == 's'),
	fid=fopen(fname,'r','ieee-be');
else
	fid=fopen(fname,'r');
end
	
ns = ns2 - ns1 + 1;
nt = nt2 - nt1 + 1;

im=zeros(nx,ny,ns,nt);
im=reshape(im,nx,ny,ns,nt);
%h=waitbar(0,'Loading images');

% Move to beginning of first time point nt1. Need to use tns here!
% fseek and ftell use bytes, need to use *2 for uint16!
if (nt1 > 1),
	fseek(fid,2*nx*ny*tns*(nt1-1),-1);
end

% added 29-08-00 to allow reading selected slices. 
% Move to beginning of first selected slice ns1:
if (tns ~= ns),
   fseek(fid,2*nx*ny*(ns1-1),0);
end


for j=1:nt,
    %waitbar(j/nt,h)
    % ftell(fid)/(2*128^2)
%    im0=fread(fid,nx*ny*ns,'int16'); % advances pointer 2*nx*ny*ns
    im0=fread(fid,nx*ny*ns,'int32'); % advances pointer 2*nx*ny*ns
    im0=reshape(im0,nx,ny,ns);
    image(im0);
    im(:,:,:,j)=im0(:,:,:);	% im0 has size(nx,ny,ns)
    % Move to beginning of next ns1:
    if (tns ~= ns),
	    fseek(fid,2*nx*ny*(tns-ns),0);
    end

end

if (size(im,3) == 1),
	im = squeeze(im);
end

fclose(fid);
%close(h);

