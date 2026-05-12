function displayTracePlot(type,tr,dat)

global DATA CUR PLT

if strcmp(type,'xcov')
  set(PLT.XsliderA(tr),'Value',CUR.Xchans(tr,1));
  ch1 = CUR.Xchans(tr,1);
  col1 = getAreaColor(CUR.areas(ch1));
  set(PLT.XeditA(tr),'String',sprintf('%d',CUR.electrodes(ch1)),'ForegroundColor',col1);

  set(PLT.XsliderB(tr),'Value',CUR.Xchans(tr,2));
  ch2 = CUR.Xchans(tr,2);
  col2 = getAreaColor(CUR.areas(ch2));
  set(PLT.XeditB(tr),'String',sprintf('%d',CUR.electrodes(ch2)),'ForegroundColor',col2);
  xcindx = getXCovIndex(CUR.Xchans(tr,1),CUR.Xchans(tr,2));  
  set(PLT.Xplot(tr),'YData',dat(:,xcindx));
  placeTitleByHandle(PLT.Xtitle(tr),0.85,0.8);
else
  ch  = CUR.chans(tr);
  col = getAreaColor(CUR.areas(ch));
  set(PLT.slider(tr),'Value',ch);
  set(PLT.edit(tr),'String',sprintf('%d',ch),'ForegroundColor',col);
  set(PLT.plot(tr),'YData',dat(:,ch));
  placeTitleByHandle(PLT.title(tr),0.85,0.8);
end
refresh(122);



%--------------
%--------------

function indx = getXCovIndex(x,y)

mx = max(x,y);
mn = min(x,y);

subsum = 0;
for i=1:mx-1
    subsum = subsum+i;
end
indx = subsum+mn;