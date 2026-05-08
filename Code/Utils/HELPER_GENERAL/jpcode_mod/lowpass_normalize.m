function [normalized] = lowpass_normalize(image, width, amplifi,devation)
%
%function lowpass_normalize
%
%call:      [normalized] =...   ; or
%           normalized = lowpass_normalize(image); oder
%           normalized = lowpass_normalize(image, width, amplifi, devation)
%
%makes a lowpass based normalization
%
%input parameter:   image       inhomogen image array
%                               optimum would be a 128x128 Matrix or a 256x256 Matrix
%                               at other sizes is filled up with zeros to the next
%                               suitable Matrix size
%                   optional:
%                   width       filter width is is parameter for the
%                               strength of the lowpass filter
%                               smal values are stading for a
%                               strong lowpass filter!
%                               width should be a multiple of 2 and in the
%                               range between ...
%                               default: ... 
%                   amplifi     maximal allowed amplification
%                               default: amplifi=5
%                   devation    parameter for the prepare smoothing of
%                               the image
%                               default= ...(30); 1000 means no smoothing and
%                               smaler value a more intensive smoothing
%                               too smal values counteract to the
%                               normalization and too big values leads to a
%                               exaggerated amplification at the borders
%
%output parameter:  normalized  normalized image array
%
%
%
%                               viel Spaß damit!
%
%                            © MB (Michael Beyerlein)
%                                 Bayo.M@gmx.de


%-Anmerkungen zum Programm selbst-

error(nargchk(1,4,nargin)) %für den Fall, dass zu wenig Input

if nargin < 2, width = 10; amplifi = 5; devation = 20; 

%preperation

%spikeelemination
smooth = spike_killer(image, devation);
%smoothing
smooth = edge_smooth(smooth, devation);


if max(size(image))<=128
    inv_image=fft2(smooth, 128, 128)
    inv_image=fftshift(inv_image);
elseif max(size(image))<=256
    inv_image=fft2(smooth, 256, 256); 
    inv_image=fftshift(inv_image);
end;
if max(size(image))>256
    error('The image is to big. The maximum matrix size is 256!');
end;
low=(size(inv_image,1)/2)-(width/2)+1;
high=(size(inv_image,1)/2)+(width/2);
low_frequency=inv_image(low:high,low:high);

gauss_window=gauss2d(width, width);
div=1/max(max(gauss_window((width/2):((width/2)+1), (width/2):((width/2)+1))));
gauss_window=gauss_window.*div;

low_frequency=low_frequency.*gauss_window;

inv_low_frequency=ifft2(low_frequency, 256, 256);
inv_low_frequency=abs(inv_low_frequency);
clear gauss_window low_frequency low high inv_image

div=1/max(max(inv_low_frequency));
characteristic=(inv_low_frequency.*div);

to_small=characteristic < (amplifi^(-1));
characteristic=(((~to_small).*characteristic)+(to_small.*(amplifi^(-1)))).^(-1);
%characteristic=characteristic.^(-1);

new_image=characteristic.*image;



zeigs(:,:,1)=image;
zeigs(:,:,2)=inv_low_frequency;
%zeigs(:,:,3)=characteristic;
%zeigs(:,:,4)=new_image;
%zeig_den_fit
zeig_den_fit(smooth);
zeig_den_fit(zeigs);
zeig_den_fit(new_image);
%zeig_den_fit(characteristic);
end;

%--------------------------------------------------------------------------
%smooth x,-x,y,-y,xy,-x-y,x-y,-xy

function [smooth] = edge_smooth (image, devation)

%x-direction

for i=1:(size(image,2)-1)
    allowed=image(:,i)-(image(:,i).*(devation/1000));
    ok=image(:,i+1)>=allowed;
    image(:,i+1)=(image(:,i+1).*ok)+(allowed.*(~ok));  
end;

%-x-direction

for i=1:(size(image,2)-1)
    ii=size(image,2)-i+1;
    allowed=image(:,ii)-(image(:,ii).*(devation/1000));
    ok=image(:,ii-1)>=allowed;
    image(:,ii-1)=(image(:,ii-1).*ok)+(allowed.*(~ok));  
end;

%y-direction

for j=1:(size(image,1)-1)
    allowed=image(j,:)-(image(j,:).*(devation/1000));
    ok=image(j+1,:)>=allowed;
    image(j+1,:)=(image(j+1,:).*ok)+(allowed.*(~ok));  
end;

%-y-direction

for j=1:(size(image,1)-1)
    jj=size(image,1)-j+1;
    allowed=image(jj,:)-(image(jj,:).*(devation/1000));
    ok=image(jj-1,:)>=allowed;
    image(jj-1,:)=(image(jj-1,:).*ok)+(allowed.*(~ok));  
end;

%smooth=image;

%x-y and -xy

%image_diag = spdiags(image);
%for j=1:(size(image_diag,1)-1)
 %   allowed=image_diag(j,:)-(image_diag(j,:).*((devation*sqrt(2))/1000));
  %  ok=image_diag(j+1,:)>=allowed;
   % image_diag(j+1,:)=(image_diag(j+1,:).*ok)+(allowed.*(~ok));  
   %end;
%for j=1:(size(image_diag,1)-1)
 %   jj=size(image_diag,1)-j+1;
  %  allowed=image_diag(jj,:)-(image_diag(jj,:).*((devation*sqrt(2))/1000));
   % ok=image_diag(jj-1,:)>=allowed;
    %image_diag(jj-1,:)=(image_diag(jj-1,:).*ok)+(allowed.*(~ok));  
    %end;
image_diag=zeros(size(image));
for j=1:(size(image,1)-1)
    d=diag(image, -j);
    for n=1:(size(d,1)-1)
        allowed=d(n)-(d(n)*((devation*sqrt(2))/1000));
        if d(n+1)<allowed;
            d(n+1)=allowed;
        end;
    end;
    for n=1:(size(d,1)-1)
        nn=size(d,1)-n+1;
        allowed=d(nn)-(d(n)*((devation*sqrt(2))/1000));
        if d(nn-1)<allowed;
            d(nn-1)=allowed;
        end;
    end;
    image_diag=diag(d,-j)+image_diag;
end;

%hier die lange Diagonale einfügen
d=diag(image);
for n=1:(size(d)-1)
        allowed=d(n)-(d(n)*((devation*sqrt(2))/1000));
        if d(n+1)<allowed;
            d(n+1)=allowed;
        end;
    end;
    for n=1:(size(d)-1)
        nn=size(d,1)-n+1;
        allowed=d(nn)-(d(n)*((devation*sqrt(2))/1000));
        if d(nn-1)<allowed;
            d(nn-1)=allowed;
        end; 
    end; 
image_diag=diag(d,0)+image_diag;

for i=1:(size(image,2)-1)
    d=diag(image, i);
    for n=1:(size(d)-1)
        allowed=d(n)-(d(n)*((devation*sqrt(2))/1000));
        if d(n+1)<allowed;
            d(n+1)=allowed;
        end;
    end;
    for n=1:(size(d)-1)
        nn=size(d,1)-n+1;
        allowed=d(nn)-(d(n)*((devation*sqrt(2))/1000));
        if d(nn-1)<allowed;
            d(nn-1)=allowed;
        end; 
    end; 
    image_diag=diag(d,i)+image_diag;
end;

%smooth = image_diag;

%xy and -x-y

%image=spdiags(spdiags(spdiags(image'))); %hängt aber von dem Verhältnis
%Zeilen zu Spalten ab
image_diag=spiegeln(image_diag);
image_diag2=zeros(size(image));
for j=1:(size(image_diag,1)-1)
    d=diag(image_diag, -j);
    for n=1:(size(d,1)-1)
        allowed=d(n)-(d(n)*((devation*sqrt(2))/1000));
        if d(n+1)<allowed;
            d(n+1)=allowed;
        end;
    end;
    for n=1:(size(d,1)-1)
        nn=size(d,1)-n+1;
        allowed=d(nn)-(d(n)*((devation*sqrt(2))/1000));
        if d(nn-1)<allowed;
            d(nn-1)=allowed;
        end;
    end;
    image_diag2=diag(d,-j)+image_diag2;
end;

%hier die lange Diagonale einfügen
d=diag(image_diag);
for n=1:(size(d)-1)
        allowed=d(n)-(d(n)*((devation*sqrt(2))/1000));
        if d(n+1)<allowed;
            d(n+1)=allowed;
        end;
    end;
    for n=1:(size(d)-1)
        nn=size(d,1)-n+1;
        allowed=d(nn)-(d(n)*((devation*sqrt(2))/1000));
        if d(nn-1)<allowed;
            d(nn-1)=allowed;
        end; 
    end; 
image_diag2=diag(d,0)+image_diag2;

for i=1:(size(image_diag,2)-1)
    d=diag(image_diag, i);
    for n=1:(size(d)-1)
        allowed=d(n)-(d(n)*((devation*sqrt(2))/1000));
        if d(n+1)<allowed;
            d(n+1)=allowed;
        end;
    end;
    for n=1:(size(d)-1)
        nn=size(d,1)-n+1;
        allowed=d(nn)-(d(n)*((devation*sqrt(2))/1000));
        if d(nn-1)<allowed;
            d(nn-1)=allowed;
        end; 
    end; 
    image_diag2=diag(d,i)+image_diag2;
end;
image_diag2=spiegeln(image_diag2);
smooth=image_diag2;

%--------------------------------------------------------------------------
%spike_killer

function [smooth] = spike_killer(image, devation)

border_up=image(1,:);
border_up=[border_up(1), border_up, border_up(size(border_up,2))];
border_down=image(size(image,1),:);
border_down=[border_down(1), border_down, border_down(size(border_down,2))];
border_left=image(:,1);
border_right=image(:,size(image,2));
bigger_image=[border_up; [border_left, image, border_right]; border_down];

first_order_neighborhood = [0 1 0; 1 0 1; 0 1 0];

for j=1:size(image,1)
    for i=1:size(image,2)
        jj=(j+1); ii=(i+1);
        biggest_neighbour=max(max(bigger_image((jj-1):(jj+1), (ii-1):(ii+1)).*first_order_neighborhood));
        if (biggest_neighbour+(image(j,i)*(devation/1000))) < (image(j,i))
            smooth(j,i)=(biggest_neighbour+(image(j,i)*(devation/1000)));
        else
            smooth(j,i)=image(j,i);
        end;
    end;
end;
    

%--------------------------------------------------------------------------
%generation of the 2d gauss window

function [mat] = gauss2d( sizeX, sizeY, sigmaX, sigmaY, meanX, meanY, cut)

if nargin < 2
	error('Error at the generation of the 2d gauss matrix');(biggest_neighbour+(bigger_image(j,i)*(devation/1000)))
	return; 
end;
if nargin < 3
	sigmaX = sizeX/5;
end;
if nargin < 4
	sigmaY = sizeY/5;
end;
if nargin < 5
	meanX = (sizeX+1)/2;
end;
if nargin < 6
	meanY = (sizeY+1)/2;
end;
if nargin < 7
	cut = 0;
end;

X = linspace(1, sizeX, sizeX)'* ones(1,sizeY);
Y = ones(1,sizeX)'   		  * linspace(1, sizeY, sizeY);
%[-sizeX/2:sizeX/2]'*ones(1,sizeX+1);
%Y = ones(1,sizeY+1)'   *[-sizeY/2:sizeY/2];

mat = exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
				+((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)))... 
            			/((sigmaX*sigmaY)^(0.5)*pi); 

if cut > 0
	maximun = max(max(mat))*cut;
	I = find(mat < maximun);
	mat(I) = 0;
end;


%--------------------------------------------------------------------------
%horizontale Spieglung
function [mirror] = spiegeln(image)

for j=1:size(image,1)
    mirror((size(image,1)-j+1),:)=image(j,:);
end;


%--------------------------------------------------------------------------
function zeig_den_fit(matrix)

%m=m/256;
%matrix=matrix/m;
slices=matrix(:,:,:,1);
for j=1:(size(matrix,3))
    slice=slices(:,:,j);
    show(:,:,1,j)=slice;
    %imagesc(slice,'CDataMapping','scaled');
    %image(slice,'CDataMapping','scaled');
end
figure;
colormap(gray(100));
h=montage(show);
%axis equal;
set(h, 'CDataMapping', 'scaled');
set(gca, 'clim', [(min(min(min(show)))) (max(max(max(show))))]);
%set(gca, 'clim', [0 (max(max(max(show))))]);

for i=1:(size(matrix,4)-1)
    
    input('nächster Parameter')
    slices=matrix(:,:,:,(i+1));
    
    for j=1:(size(matrix,3))
        slice=slices(:,:,j);
        show(:,:,1,j)=slice;
    end
    figure
    colormap(gray(100));
    h=montage(show);
    %axis equal;
    set(h, 'CDataMapping', 'scaled');
    set(gca, 'clim', [(min(min(min(show)))) (max(max(max(show))))]);
end;
