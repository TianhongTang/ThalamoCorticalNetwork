function displayXCov(dat, electrodes, areas, clim, plottitle)

tmp   = size(dat,3);
nchan = elecTri2Sqr(tmp);

if ~exist('electrodes') | isempty(electrodes)
  electrodes = [1:nchan];
end
if  ~exist('areas') | isempty(areas)
  areas = -1*ones(1,nchan);
end

if  ~exist('clim') | isempty(clim)
  clim = [-0.5 1.0]
end



ymin	 = -0.5;
ymax     = 1.0;
tmin	 = -128;
tmax     = 128;
len      = size(dat,1);
hlen     = (len-1)/2;
midpt    = hlen + 1;
times    = ([1:len]-midpt);
totplots = 0;

for j=1:nchan
   jch    = j;	
   elec_j = electrodes(jch);
   area_j = areas(jch);
   col_j  = getUniqueColor(elec_j,area_j);
   axlabel (0.01, 1-(nchan*0.005)-(j)/(1.2*nchan),...
	   sprintf('EL %2d (V%d)',elec_j, area_j),col_j,6,'bold',2.0);

   for i=1:nchan
    ich 	= i;
    elec_i 	= electrodes(ich);
    area_i 	= areas(ich);
    col_i  	= getUniqueColor(elec_i,area_i);
    if j == 1	
      axlabel (0.1-(0.5/nchan)+i/(1.25*nchan), 0.01,...
	       sprintf('EL %2d',elec_i),col_i,6,'bold',2.0);
    end
    if i <= j	
      idx = (j-1)*nchan+i;
      totplots = totplots + 1;
      subplot(nchan,nchan,idx);
      chidx = getTriangularIndex(jch,ich);
      imagesc(squeeze(dat(:,:,chidx)),clim);
      axis off
    end	
  end  
end

label (0.5, 0.95, plottitle);


%--------------------------------------
function idx = getTriangularIndex(i,j)
idx = sum(1:i-1)+j;

%------------------------------------------
function nchan = elecTri2Sqr(n)
vals = [1 3 6 10 15 21 28 36 45 55 66 78 91 105 120];
nchan   = find(vals == n);
