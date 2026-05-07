function x = runningaverage(x,width)
% y = runningaverage(x,width)
% calculates running average of each row using a sliding window
%  of width 'width'
% should be linear time, much faster than using convn
% twice as fast and takes less space than runningaveragenan, 
%  because doesn't try to ignore nans

assert(mod(width,2)==1,'width must be odd');

if size(x,2) < floor(width/2)
    x = repmat(mean(x,2),1,size(x,2));
    return;
end
x = [zeros(size(x,1),ceil(width/2)) x zeros(size(x,1),ceil(width/2))];
x = cumsum(x,2);
x = x(:,(width+1):(end-1)) - x(:,1:(end-width-1));

% normalize edges to account for missing data
for w = 1:floor(width/2)
    normfactor = (width/(width-floor(width/2)+w-1));
    x(:,w) = x(:,w).*normfactor;
    x(:,end-w+1) = x(:,end-w+1).*normfactor;
end

x = x ./ width;
