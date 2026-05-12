function displayXCov(dat, electrodes, areas, frm)

tmp   = size(dat,3);
nchan = elecTri2Sqr(tmp);

if ~exist('electrodes') | isempty(electrodes)
  electrodes = [1:nchan];
end
if  ~exist('areas') | isempty(areas)
  areas = -1*ones(1,nchan);
end

if  ~exist('frm') | isempty(frm)
  frm = 1;
end


ymin	= -0.5;
ymax   = 1.0;
tmin	= -128;
tmax   = 128;
len   = size(dat,1);
hlen  = (len-1)/2;
midpt = hlen + 1;
times   = ([1:len]-midpt);

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
      subplot(nchan,nchan,idx);
      if i == j
	frame_i = 'k';
	frame_j = 'k';
	lwid	  = 2;
	else
	  frame_i = col_i;
	  frame_j = col_j;	
	  lwid 	  = 1;
      end
      rectangle('Position',...
		[tmin ymin (tmax-tmin)/2 ymax-ymin],...
		'EdgeColor',frame_i,...
		'LineWidth',lwid,'FaceColor',col_i.*...
		[0.3 0.3 0.3]+get(gcf,'Color')-[0.2 0.2 0.2]);
      hold on
      rectangle('Position',...
		[tmin+(tmax-tmin)/2 ymin (tmax-tmin)/2 ymax-ymin],...
		'EdgeColor',frame_j, 'LineWidth',lwid,'FaceColor',col_j.*...
		[0.3 0.3 0.3]+get(gcf,'Color')-[0.2 0.2 0.2]);

      chidx = getTriangularIndex(jch,ich);
      tmpwv = mean(dat(:,frm,chidx),2);
      plot(times,tmpwv);
      line([0 0],[ymin ymax],'Color','k');
      line([tmin tmax],[0 0],'Color','k');
      set(gca,'XLim',[tmin tmax],'YLim',[ymin ymax]);
      axis off
    end	
  end  
end



%--------------------------------------
function idx = getTriangularIndex(i,j)
idx = sum(1:i-1)+j;

%------------------------------------------
function nchan = elecTri2Sqr(n)
vals = [1 3 6 10 15 21 28 36 45 55 66 78 91 105 120];
nchan   = find(vals == n);
