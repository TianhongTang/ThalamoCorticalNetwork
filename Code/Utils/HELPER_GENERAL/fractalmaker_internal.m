function [IMG]= fractalmaker_internal

param = struct();
rand('twister',sum(100*clock));
kuikomisize=2;
dekoboko=0.8;
param.hires_image_resolution_dpi = 150;
param.lores_image_number_of_pixels = 150
numofsuperimp=4;
datamat=[];
for sup=[1:numofsuperimp]
    numofedge= 4+round(rand*6);
    edgesize=numofsuperimp-sup+rand;
    numofrecursion= 2+round(rand*3);
    fracx=edgesize*cos(2*pi*[1:numofedge]/numofedge);
    fracy=edgesize*sin(2*pi*[1:numofedge]/numofedge);
    GAmat=[];
    for k=[1:1:numofrecursion]
        mx=(fracx([2:end 1])+fracx)/2;
        my=(fracy([2:end 1])+fracy)/2;
        dx= fracx([2:end 1])-fracx;
        dy= fracy([2:end 1])-fracy;
        theta=atan(dy./dx);
        theta(find(dx<0))=theta(find(dx<0))+pi;
        GA=kuikomisize*(rand-dekoboko);
        fracx2=mx+GA*sin(theta);
        fracy2=my-GA*cos(theta);
        fracx=[fracx;fracx2];
        fracx=fracx(:);
        fracx=fracx';
        fracy=[fracy;fracy2];
        fracy=fracy(:);
        fracy=fracy';
        GAmat=[GAmat;GA];
    end
    col=rand(1,3);
    
%uncomment for isoluminance
%     col(1)=0;
%     [tz]=randperm(2);
%     nega=[-1 1]; nega=nega(tz(1));
%     col=col*nega;
%     [tz]=randperm(2);
%     nega=[-1 1]; nega=nega(tz(1));
%     col=col*nega;
%     col=dkl2rgb(col')';
    
    if sup==1
    h=figure('Visible','Off','PaperPositionMode', 'auto')
    end
    fill(fracx([1:end 1]),fracy([1:end 1]),col)
    hold on
    plot(fracx([1:end 1]),fracy([1:end 1]),'Color',col)
    axis square
    axis([-numofsuperimp-1 numofsuperimp+1 -numofsuperimp-1 numofsuperimp+1])
    hold on
    

end

axis square
axis([-numofsuperimp-1 numofsuperimp+1 -numofsuperimp-1 numofsuperimp+1])
axis([-5,5 -5,5])
axis off
set(gcf,'color','k','inverthardcopy','off');
set(gcf,'menubar','none');
set(gca,'Units','pixels');
posa=get(gca,'Pos');
set(gca,'Units','normalized');
posf=get(gcf,'Pos');
set(gcf,'Pos',[posf(1) posf(2) posa(3) posa(3)]);
set(gcf,'paperpositionmode','auto');
%IMG=frame2im(getframe(h)); %this is way faster than print in high res..
IMG = print(h,'-RGBImage',['-r' num2str(param.hires_image_resolution_dpi)]);
clf;
close(h); 
IMG = trim_image(IMG,'symmetrical','center','leaveborder',1,'resize',{param.lores_image_number_of_pixels.*[1 1],'lanczos3'});





