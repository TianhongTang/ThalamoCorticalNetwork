function [h] = plot_errorpatch(x,y,lo,hi,varargin)
% [h] = plot_errorpatch(x,y,lo,hi,varargin)
% draw a line on top of a patch; 
%  the patch extends between the lower and upper error bounds
%
% if varargin gets two cell arrays as arguments, uses the first
%  as the style for the main plot line, and the second as
%  the style for the patch.
%
%
% examples:
%  
%  x = 1:100;
%  y = randn(1,100);
%  yse = 1 + abs(randn(1,100));
%
%
%  % plot the data as a red line with small dots, 
%  % and the error bounds as a gray patch with black edges.
%  % use the vector yse to provide the standard error
%  % for each data point:
%  plot_errorpatch(x,y,yse,{'r.-'},{'facecolor',[.5 .5 .5],'edgecolor','k'});
%
%  % the same as above, except separately specifying the lower bound (y-yse)
%  % and the upper bound (y+yse) for each data point
%  % by manually specifying the lower and upper bounds, you can plot
%  % error bounds even when they are asymmetric (e.g. for bootstrap 
%  % confidence intervals)
%  plot_errorpatch(x,y,y-yse,y+yse,{'r.-'},{'facecolor',[.5 .5 .5],'edgecolor','k'});

if ~isvector(x) || ~isvector(y)
    error('inputs must be vectors');
end;

if nargin < 3
    error('too few args');
elseif nargin == 3
    hi = y + lo;
    lo = y - lo;
elseif nargin >= 4
    if ischar(hi) || iscell(hi)
        varargin = [{hi} varargin];
        hi = y + lo;
        lo = y - lo;
    end;
end;

x = torowvector(x);
y = torowvector(y);
lo = torowvector(lo);
hi = torowvector(hi);

use_separate_styles = length(varargin) == 2 && iscell(varargin{1}) && iscell(varargin{2});
if use_separate_styles
    arg1 = varargin{1};
    arg2 = varargin{2};
else
    arg1 = varargin;
    arg2 = varargin;
end;

% sadly, we have to specify a color to 'patch' before listing variable
%  args. So if a color isn't already being specified in the variable args,
%  then we have to choose a default color 
% (we choose the line's color, made slightly lighter)
patch_color_not_specified = ~use_separate_styles &&  ( (length(arg2) < 1) || (~hasparam(arg2,'FaceColor') && ~iscolorspec(arg2{1})) );
patch_color_specified_immediately = ~patch_color_not_specified && iscolorspec(arg2{1});
default_patch_color = [.7 .7 .7];

was_hold = ishold;
hold on;

% plot a separate patch for each run of non-NaN standard errors
curnan = isnan(lo) | isnan(hi);
if any(curnan)
    h = [];
    curstart = find(~curnan & prev(curnan,1,true));
    curend = find(~curnan & next(curnan,1,true));
    for i =1:numel(curstart)
        h = [h ; plot_errorpatch( ...
            x(curstart:curend), ...
            y(curstart:curend), ...
            lo(curstart:curend), ...
            hi(curstart:curend), ...
            arg1,arg2)];
    end;
else
    curstart = 1;
    curend = numel(x);
end;

if patch_color_not_specified
    % make sure that Patch can handle the color spec
    for a = 1:length(arg2)
        if ischar(arg2{a}) && size(arg2{a},1)==1 && ~isempty(strmatch('co',lower(arg2{a})))
            arg2{a} = 'FaceColor';
        end;
    end;
    

    % plotpatch(i) = plot the i-th patch
    plotpatch2 = @(x,lo,hi) patch([x x(end:-1:1)],[lo hi(end:-1:1)],default_patch_color,arg2{:});
    plotpatch = @(i) plotpatch2(x(curstart(i):curend(i)),lo(curstart(i):curend(i)),hi(curstart(i):curend(i)));
    h2 = [];
    for i = 1:numel(curstart)
        h2 = [h2 ; plotpatch(i)];
    end;
    h1 = plot(x,y,arg1{:});
    
    linecolor = colorspec_to_rgb(get(h1,'Color'));
    patchcolor = interpcolor(linecolor,[1 1 1],.5);
    set(h2,'FaceColor',patchcolor);
elseif patch_color_specified_immediately
    % plotpatch(i) = plot the i-th patch
    plotpatch2 = @(x,lo,hi) patch([x x(end:-1:1)],[lo hi(end:-1:1)],arg2{:});
    plotpatch = @(i) plotpatch2(x(curstart(i):curend(i)),lo(curstart(i):curend(i)),hi(curstart(i):curend(i)));
    h2 = [];
    for i = 1:numel(curstart)
        h2 = [h2 ; plotpatch(i)];
    end;
    h1 = plot(x,y,arg1{:});
else
    % plotpatch(i) = plot the i-th patch
    plotpatch2 = @(x,lo,hi) patch([x x(end:-1:1)],[lo hi(end:-1:1)],default_patch_color,arg2{:});
    plotpatch = @(i) plotpatch2(x(curstart(i):curend(i)),lo(curstart(i):curend(i)),hi(curstart(i):curend(i)));
    h2 = [];
    for i = 1:numel(curstart)
        h2 = [h2 ; plotpatch(i)];
    end;
    h1 = plot(x,y,arg1{:});
end;

% default: no edge color on patch, unless specified by user
if ~hasparam(arg2,'EdgeColor')
    set(h2,'EdgeColor','none');
end;

h = [h1 ; h2];

if ~was_hold
    hold off;
end;
