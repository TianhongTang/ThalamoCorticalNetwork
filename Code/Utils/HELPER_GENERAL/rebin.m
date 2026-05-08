function [y] = rebin(x,binwidth,combine_fxn,shrink,truncate)
%rebin(x,binwidth,combine_fxn,shrink,truncate)
%
% bins data by combining 'binwidth' adjacent data points
% into a single new data point, using the combine_fxn (default=@mean)
% if shrink is true, just return these new data points
% if shrink is false (the default), return a data vector
%  of the same size as x, with each
%  data point replaced by its new data point
% if truncate is true or 'end', then if x doesn't divide evenly into
%  bins of width 'binwidth', then remove the remainder of x.
% if truncate is 'start', truncate the starting entries of x
%  instead of the ending entries of x.
% if truncate is 'startnan' or 'endnan', instead of removing the
%  entries, set them to nan
%
% if applied to a matrix, is applied separately to each row
%
% Update 2018-05-27: added much faster computation for the special case
% where the following conditions are met:
%  combine_fxn = @sum or @mean
%  shrink = true,
%  truncate = false
% Note: this does not not actually invoke a call to "sum" or "mean",
% instead it just calls "cumsum".
%
%examples
% rebin([1 2 3 4],2,@mean) == [1.5 1.5 3.5 3.5]
% rebin([1 2 3 4],2,@mean,true) == [1.5 3.5]
% rebin([1 2 3 4],2,@sum) == [3 3 7 7]
% rebin([1 2 3 4],2,@sum,true) == [3 7]
% rebin([1 2 3 4 5],2,@sum,true,true) == [3 7]
% rebin([1 2 3 4; 5 6 7 8],2,@mean,true) == [1.5 3.5 ; 5.5 7.5]

if nargin < 1 error('need at least 1 argument'); end
if nargin < 2 binwidth = 1; end
if nargin < 3 combine_fxn = @mean; end
if nargin < 4 shrink = false; end
if nargin < 5 truncate = false; end
if nargin > 5 error('can take at most 5 arguments'); end

if isequal(truncate,'start') || isequal(truncate,'end') || isequal(truncate,'startnan') || isequal(truncate,'endnan')
    do_truncate = true;
elseif isequal(truncate,true)
    do_truncate = true;
    truncate = 'end';
else
    do_truncate = false;
    truncate = false;
end

if ~do_truncate && mod(size(x,2),binwidth) ~= 0
    error(['data does not divide evenly into bins of width ' num2str(binwidth)]);
end


% special case for sum/mean or nansum/nanmean while shrinking, 
% which can be executed very fast using "cumsum" instead of sum/mean 
% functions
if shrink && ~do_truncate 
    
    % nansum or nanmean?
    if isequal(combine_fxn,@nansum) || isequal(combine_fxn,@nanmean)
        x_is_nan = isnan(x);
        total_num_nan = sum(x_is_nan(:));
        if total_num_nan == 0
            % if there are no NaNs in the data, we can use the more
            % efficient algorithm for doing the simple sum or mean.
            if isequal(combine_fxn,@nansum)
                combine_fxn = @sum;
            else
                combine_fxn = @mean;
            end
        else
            % use the algorithm for handling NaNs:
            % set all NaNs to zero, do the calculation as usual, 
            % and then later correct for the missing entries
            x_with_nans_zeroed = x;
            x_with_nans_zeroed(x_is_nan) = 0;

            % sum of x (with NaNs ignored) in each bin
            csum = cumsum(x_with_nans_zeroed,2);
            y = csum(:,binwidth:binwidth:end) ...
                - csum(:,1:binwidth:end) ...
                + x_with_nans_zeroed(:,1:binwidth:end);

            csum_nan = cumsum(x_is_nan,2);
            numnan = csum_nan(:,binwidth:binwidth:end) ...
                - csum_nan(:,1:binwidth:end) ...
                + x_is_nan(:,1:binwidth:end);

            % if is sum, set to NaN for bins where all data points are NaN
            % if is mean, get the mean by dividing by the number of non-NaN 
            %  data points per bin (which again will cause it to be set to NaN
            %  if all data points in the bin are NaN, since 0/0 = NaN).
            if isequal(combine_fxn,@nansum)
                y(numnan == binwidth) = nan;
            else
                y = y ./ (binwidth - numnan);
            end

            return;
        end
    end
        
    % sum or mean?
    if isequal(combine_fxn,@sum) || isequal(combine_fxn,@mean)
        % calculate sum in each bin
        csum = cumsum(x,2);
        y = csum(:,binwidth:binwidth:end) ...
             - csum(:,1:binwidth:end) ...
             + x(:,1:binwidth:end);

        % if needed, get the mean by dividing by the number of data points per
        % bin
        if isequal(combine_fxn,@mean)
            y = y ./ binwidth;
        end

        return;
        
    else

    end
end


% if has multiple rows...
if length(size(x)) > 2 || size(x,1) > 1 && size(x,2) > 1
    % rebin separately on each row
    y_init = rebin(x(1,:),binwidth,combine_fxn,shrink,truncate);
    y = nans(size(x,1),size(y_init,2));
    y(1,:) = y_init;
    for r = 2:size(x,1)
        y(r,:) = rebin(x(r,:),binwidth,combine_fxn,shrink,truncate);
    end
    return;
end

nbins = floor(length(x) / binwidth);

if ~do_truncate || isequal(truncate,'end') || isequal(truncate,'endnan')
    binstart = 1:binwidth:length(x);
    if binstart(end)+binwidth-1 > length(x)
        binstart = binstart(1:(end-1));
    end
    binend = binstart + binwidth - 1;
elseif isequal(truncate,'start') || isequal(truncate,'startnan')
    binend = length(x):-binwidth:1;
    if binend(end) < binwidth
        binend = binend(1:(end-1));
    end
    binend = binend(end:-1:1);
    binstart = binend - binwidth + 1;
else
    error('unknown truncation mode!');
end

if shrink
    if size(x,1) > size(x,2)
        y = zeros(nbins,1);
    else
        y = zeros(1,nbins);
    end
    for b = 1:nbins
        y(b) = combine_fxn(x(binstart(b):binend(b)));
    end
else
    y = nans(size(x));
    for b = 1:nbins
        y(binstart(b):binend(b)) = combine_fxn(x(binstart(b):binend(b)));
    end
    
    if isequal(truncate,'end') || isequal(truncate,'start')
        y = y(binstart(1):binend(end));
    end
end
