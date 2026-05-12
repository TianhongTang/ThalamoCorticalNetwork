function hout = offsetplot(offset, varargin)
% OFFSETPLOT  Plot a series of lines with offsets
% Usage:
%    offsetplot(offset,Y)
%    offsetplot(offset,X,Yl,...)
%    offsetplot(offset,X,Yl,LineSpec,...)
%    offsetplot(...,'PropertyName',PropertyValue,...)
%    offsetplot(offset,axes_handle,...)
%    offsetplot([],...)
%    h = offsetplot(...)
%
% OFFSETPLOT plots a series of curves using PLOT, but instead of overlaying
% them on top of one another, it applies an offset to each curve. The
% resulting plot, with several "parallel" curves, is a standard way of
% showing, for example, EEG traces.
%
% OFFSETPLOT can be used just as PLOT is used, only adding the offset
% argument to the begining of the argument list. All arguments supplied
% after the offset argument are passed to PLOT without modification.
%
% The argument "offset" can be a scalar, in which case the first curve in
% the series is plotted with no offset, the second curve is moved by an
% amount -offset, the third is moved by -2*offset, and so on. This produces
% a series of curves with the first curve at the top of the graph, and
% subsequent curves appearing lower on the graph. For the opposite
% behavior, use a negative value for "offset". "offset" can also be a
% vector of length N, where N is the number of curves to be created. In
% this case, the offsets are applied element-wise from the vector.
%
% "offset" can also take any of the following special values:
%     [], NaN, Inf, -Inf
% In this case, the offset is automatically calculated to produce the
% minimum value necessary to guarantee that adjacent curves do not overlap,
% or four times the maximum standard deviation of all the curves, whichever
% is larger. In the case of -Inf, the graph is generated with the first
% curve lowest, with subsequent curves ascending; otherwise, the first
% curve is highest on the graph, with subsequent curves descending.
%
% Example: Using default offset behavior
%    x = 0:.1:10;
%    Y = randn(101,10);
%    offsetplot([],x,Y,'--k')
% Example: Specify a uniform offset
%    t   = (0:.01:1)';
%    phi = 0:pi/5:pi;
%    f   = 2;
%    T   = repmat(t,size(phi));
%    PHI = repmat(phi,size(t));
%    Y   = sin(2*pi*T*f + PHI);
%    offsetplot(1,t,Y,'r','LineWidth',2)
% Example: Specify a vector offset
%    x = 0:.1:10;
%    Y = randn(101,5);
%    off = [0 10 -10 20 -20];
%    offsetplot(off,x,Y);

% Plot the lines as requested
h = plot(varargin{:});
if nargout > 0
    hout = h; % Return an output only if one is requested
end

% To improve portability, we will pretend we know nothing at all about the
% input interface to PLOT. This means that this code should automatically
% support any new additions to PLOT that MATLAB adds in the future. We will
% get all the information we need to know from the line handles returned by
% PLOT.
nlines = numel(h);
if nlines<2
    return;
end
X = get(h,{'XData'});
if ~congruent(X)
    % At the moment, we will only support all lines having the same XData.
    error('offsetplot:nonuniformX', ...
        'XData must be the same for all lines.');
end
Y = get(h,{'YData'});
% Convert the cell array into a matrix (this is a slightly surer bet than
% using cell2mat, because we don't necessarily know whether any given YData
% is a row or column vector):
Y = cellfun(@(x)(x(:)),Y,'UniformOutput',false);
Yy = [Y{:}];
if isempty(offset) ...
        || (isscalar(offset) && (offset==Inf||offset==-Inf||isnan(offset)))
    % Use the minimum offset that guarantees no overlap, or four times the
    % maximum standard deviation, whichever is larger.
    overlap = diff(Yy,1,2);
    max_overlap = max(overlap(:));
    std_range = 4*max(std(Yy));
    if offset == -Inf
        coef = -1; % Move up instead of down
    else
        coef = 1;
    end
    offset = coef * max([max_overlap std_range]);
end

% Generate the offset vector:
offset_mult = 0:-1:(1-nlines);
offset = offset(:)'.*offset_mult; % This is the line that will produce an
                                  % error if the user supplies a vector
                                  % with the wrong number of elements.
axh = get(h(1),'Parent'); % Get the axes the lines were plotted onto
% Calculate the new axes limits:
ylimit = get(axh,'YLim');
offsetlimit = [min(offset) max(offset)];
new_ylimit = ylimit + offsetlimit;
% Calculate the new tick positions. We want to put one at every zero-value
% to show the midpoints, and one at the top and bottom of the graph to show
% the range in the original units.
ytick = sort(offset(end:-1:1));
ytick = [ylimit(1)+ytick(1) ytick ylimit(2)+ytick(end)];
ytick_labs = [ylimit(1) zeros(size(offset)) ylimit(2)];

% Move the curves.
for k=1:nlines
    set(h(k),'YData',Y{k}+offset(k));
end
set(axh,'YLim',new_ylimit,'YTick',ytick,'YTickLabel',ytick_labs);

function tf = congruent(C)
% CONGRUENT Test for congruency of all cell array elements
% Returns true if all elements in a cell array are exactly identical, false
% otherwise.
ns = cellfun(@numel,C);
tf = all(ns==ns(1));
if ~tf
    return;
end
tf = all(cellfun(@(x)(all(x==C{1})),C));
