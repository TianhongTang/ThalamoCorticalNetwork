function fullplx = iplx_mimichelper(filename,closeFlag)

persistent ifullplx;

if nargin < 2
    closeFlag = false;
end

if isempty(filename)
    if closeFlag
        ifullplx = [];
        fullplx = [];
        return;
    else
        [filename,pathname] = uigetfile('*.plx', 'Choose a PLX file');
        filename = fullfile(pathname, filename);
    end
end

if ~isempty(ifullplx) && strcmp(ifullplx.Filename,filename)
    if closeFlag
        ifullplx = [];
    end
else
    ifullplx = iplx_fullfile(filename);
end

fullplx = ifullplx;