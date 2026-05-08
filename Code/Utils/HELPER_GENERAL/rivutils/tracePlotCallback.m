function tracePlotCallback(tag)

global CUR

type = 'normal';

if ~exist('tag') | isempty(tag)
   obj  = gcbo;
   tag  = get(obj,'Tag');
   val  = round(get(obj,'Value'));
   dat   = get(obj,'UserData');
   str  = get(obj,'String');
end
indx = findstr(tag,'_');
tr   = str2num(tag(indx+1));
token = tag(1:indx-1);
if strcmp(token,'Slider')
    CUR.chans(tr) = round(val);
elseif strcmp(token,'SliderA')
    CUR.Xchans(tr,1) = round(val);
    type = 'xcov';
elseif strcmp(token,'SliderB')
    CUR.Xchans(tr,2) = round(val);
    type = 'xcov';
elseif strcmp(token,'Edit')
    CUR.chans(tr) = round(str2num(str))
elseif strcmp(token,'EditA')
    CUR.Xchans(tr,1) = round(str2num(str))
    type = 'xcov';
elseif strcmp(token,'EditB')
    CUR.Xchans(tr,2) = round(str2num(str));
    type = 'xcov';
end

displayTracePlot(type,tr,dat);
