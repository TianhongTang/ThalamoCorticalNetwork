function varargout = errorfield(varargin)
% ERRORFIELD Add a colored region representing uncertainty to a plot
%    errorfield(Y,E)
%    errorfield(X,Y,E)
%    errorfield(X,Y,L,U)
%    errorfield(haxes,...)
%    errorfield(hline,E)
%    errorfield(hline,L,U)
%    errorfield(...,'PropertyName',PropertyValue,...)
%    hpatch = errorfield(...)
%  The behavior of ERRORFIELD is similar to that of ERRORBAR, but instead
%  of drawing error bars at each point in Y, it draws a patch object
%  to produce a colored field around a line representing the uncertainty.
%  Unlike ERRORBAR, ERRORFIELD does not draw the line itself; this is to
%  allow flexibility in case plotting the line is undesired.
%
%  ERRORFIELD differs from ERRORBAR also in allowing cell arrays to be
%  passed as arguments for X, Y, E, L, and U. In these cases, the cell
%  arrays should all have the same length. Additionally, the k'th element
%  of each cell array should be a vector such that any supplied X{k}, Y{k},
%  E{k}, L{k}, or U{k} should have the same length. This allows errorfields
%  around lines of differing numbers of points.
%  
%  The usage
%    errorfield(...,'PropertyName',PropertyValue,...)
%  allows additional properties to be passed to the patch object. The
%  additional 'PropertyName' option 'Color' is also defined, which takes a
%  ColorSpec argument to set the errorfield's color. If multiple
%  errorfields are being defined, each PropertyValue may optionally be a
%  cell vector of values, one for each errorfield.
%  
%  The usage
%    errorfield(haxes,...)
%  plots the errorfield in the axes whose handle is h.
%  
%  The usages
%    errorfield(hline,E,...)
%    errorfield(hline,L,U,...)
%  plot the errorfield around the line(s) given by the handle(s) hline.
%  In this case, the length of E or of L and U must be the same as the
%  length of the line's XData and YData properties. See the examples below.
%  
%  Example: Plotting a line with an errorfield
%    plot(x,y,'b')
%    errorfield(x,y,e,'Color','b')
%  Example: Plotting a line with an errorfield, with automatic coloration
%    hline = plot(x,y)
%    errorfield(hline,edown,eup);
%  Example: Showing off the best features of ERRORFIELD all at once
%    errorfield(plot(x1,y1,x2,y2),{e1,e2},'FaceAlpha',{.75,.25})
%
%  See also ERRORBAR, PLOT, PATCH.

% Geoffrey Adams 8/16/07

if nargin < 1
    error('errorfield:inputs', 'Insufficient number of input arguments.');
end

% Set defaults
plotAxes = gca; % Axes to plot on
numfields = 1;  % Number of errorfields we are plotting
linehandles = false; % Whether or not the first argument is line handles
color = {}; % Errorfield ColorSpec values.

% Find character inputs -- they should be the PropertyName arguments. If
% there are none, then there are no properties passed. Otherwise,
% everything before the first one is either data or handle.
chars = find(cellfun(@ischar,varargin));
if isempty(chars)
    props = {};
    datas = varargin;
else
    firstProp = chars(1);
    props = varargin(firstProp:end);
    datas = varargin(1:firstProp-1);
end

% Check for handle inputs
if all(ishandle(datas{1}))
    % The first argument is a handle, so it's not data.
    h = datas{1};
    datas = datas(2:end);
    % We will be plotting as many errorfields as there are handles:
    numfields = numel(h);
    hType = get(h,'Type');
    if all(strcmp('axes', hType))
        % The handles passed were axes handles, so we're going to plot in
        % those axes.
        plotAxes = h;
    elseif all(strcmp('line', hType));
        % The handles passed were line handles, so we should pick up the
        % data and color information from the lines.
        linehandles = true;
        plotAxes = get(h,'Parent');
        if numfields > 1
            % If there are more than one line handle, we need to massage
            % things a bit to get it into the right format:
            plotAxes = cell2mat(plotAxes);
            %XCell = get(h,'XData');
            %YCell = get(h,'YData');
            X = get(h,'XData');
            Y = get(h,'YData');
            color = get(h,'Color');
            %newshape = repmat({[numfields,1]},size(XCell));
            %XCell = cellfun(@reshape,XCell,newshape);
            %YCell = cellfun(@reshape,YCell,newshape);
            %X = [XCell{:}];
            %Y = [YCell{:}];
        else
            % No such problem with only one handle:
            X = {get(h,'XData')};
            Y = {get(h,'YData')};
            color = {get(h,'Color')};
        end
%         if isvector(X)
%             X = X(:);
%         end
%         if isvector(Y)
%             Y = Y(:);
%         end
    else
        error('errorfield:handles', ['Leading handle arguments must ', ...
            'be all axes handles or all line handles.']);
    end
end

if ~linehandles
    % If the a linehandle wasn't passed, then Y and possibly X arguments
    % must have been:
    switch numel(datas)
        case {0,1}
            error('errorfield:inputs', ['Insufficient number of ' ...
                'numeric inputs.']);
        case 2
            % Calling syntax was something of the form
            % errorbar(Y,...)
            % (or optionally with an haxes argument)
            Ymat = datas{1};
            if iscell(Ymat)
                Y = Ymat;
            else
                if ~isvector(Ymat)
                    Y = mat2cell(Ymat, size(Ymat,1), ones(1,size(Ymat,2)));
                else
                    Y = {Ymat(:)};
                end
            end
            if numfields > 1 && numfields ~= size(Y,2)
                error('errorfield:dimMismatch', ['Number of Y plots ' ...
                    'must match number of supplied axes handles.']);
            end
            X = cellfun(@defaultvector,Y);
            datas = datas(2);
        case {3,4}
            % Calling syntax was something of the form
            % errorbar(X,Y,...)
            % (or optionally with an haxes argument)
            Xmat = datas{1};
            Ymat = datas{2};
            if iscell(Xmat)
                X = Xmat;
            else
                if ~isvector(Xmat)
                    X = mat2cell(Xmat, size(Xmat,1), ones(1,size(Xmat,2)));
                else
                    X = {Xmat};
                end
            end
            if iscell(Ymat)
                Y = Ymat;
            else
                if ~isvector(Ymat)
                    Y = mat2cell(Ymat, size(Ymat,1), ones(1,size(Ymat,2)));
                else
                    Y = {Ymat};
                end
            end
%             if isvector(X)
%                 X = X(:);
%             end
%             if isvector(Y)
%                 Y = Y(:);
%             end
            if ~(numel(X)==numel(Y) && all(cellfun(@samesize,X,Y)))
                error('errorfield:dimMismatch', ['X and Y must be the ' ...
                    'same size.']);
            elseif numfields > 1 && numfields ~= numel(Y)
                error('errorfield:dimMismatch', ['Number of Y plots ' ...
                    'must match number of supplied axes handles.']);
            end
            datas = datas(3:end);
        otherwise
            error('errorfield:inputs', 'Too many leading numeric inputs.');
    end
end

X = cellfun(@makecolumnvector,X,'UniformOutput',false);
Y = cellfun(@makecolumnvector,Y,'UniformOutput',false);

numfields = numel(Y);
if numel(plotAxes) < numfields
    plotAxes = repmat(plotAxes,[numfields,1]);
end

switch numel(datas)
    case 1
        Lmat = datas{1};
        Umat = datas{1};
    case 2
        Lmat = datas{1};
        Umat = datas{2};
end
if iscell(Lmat)
    L = Lmat;
else
    if ~isvector(Lmat)
        L = mat2cell(Lmat, size(Lmat,1), ones(1,size(Lmat,2)));
    else
        L = {Lmat(:)};
    end
end
if iscell(Umat)
    U = Umat;
else
    if ~isvector(Umat)
        U = mat2cell(Umat, size(Umat,1), ones(1,size(Umat,2)));
    else
        U = {Umat(:)};
    end
end

oldprops = props;
for k=1:2:numel(oldprops)
    if strcmpi('color',oldprops{k})
        color = oldprops{k+1};
        if ~iscell(color)
            color = {color};
        end
        props = props([1:k-1, k+2:end]);
    end
end

vx = cell(numfields,1);
vy = vx;
for k=1:numfields
    vx{k} = [X{k}; X{k}(end:-1:1,:)];
    vy{k} = [(Y{k}-L{k}); (Y{k}(end:-1:1,:)+U{k}(end:-1:1,:))];
end

if isempty(color)
    color = repmat({'k'},[numfields,1]);
elseif isscalar(color)
    color = repmat(color,[numfields,1]);
end

thisprops = props;
if ~isempty(props)
    cellprops = cellfun(@iscell,props);
    scalarprops = cellfun(@isscalar,props);
    for k=find(cellprops&scalarprops)
        thisprops{k} = thisprops{k}{1};
    end
    cellprops = find(cellprops&~scalarprops);
else
    cellprops = [];
end

hpatch = zeros(numfields,1);

for k=1:numfields
    thisvx = vx{k};
    thisvy = vy{k};
    thiscolor = color{k};
    thisaxes = plotAxes(k);
    for m=cellprops
        thisprops{m} = props{m}{k};
    end
    hpatch(k) = patch(thisvx, thisvy, thiscolor, ...
        'Parent', thisaxes, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
        thisprops{:});
end

switch nargout
    case 0
        % Do nothing
    case 1
        varargout{1} = hpatch;
    otherwise
        error('errorfield:outputs', 'Too many output arguments.');
end
    
function tf = samesize(A,B)
if isvector(A) && isvector(B)
    tf = numel(A)==numel(B);
else
    tf = all(size(A)==size(B));
end

function x = defaultvector(y)
x = 1:numel(y);

function a = makecolumnvector(a)
if isvector(a)
    a = a(:);
end