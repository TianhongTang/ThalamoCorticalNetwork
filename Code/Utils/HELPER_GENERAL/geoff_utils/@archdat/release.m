function ob = release(ob)
% ARCHDAT/RELEASE  Release archdat data from working memory
% Usage:
%   ad = release(ad)
% archdat objects refer to data on disk, but in some cases it may be
% desirable to hold the data in memory for a while. This "readies" the data
% for repeated usage, and can be done with the READY method. After being
% readied, the archdat can release the data from memory with the RELEASE
% method. An example:
%    ad = ready(ad);
%    for k=1:n
%        ...
%        someFunction(ad);
%        ...
%    end
%    ad = release(ad);
% In this example, someFunction takes an archdat object as an argument.
% Inside the function, it may call methods that require disk access, such
% as GET, DATASIZE, DATACLASS, or INFO. It would be inefficient to access
% the disk repeatedly every time the function is called, so (provided
% enough memory is available) it is better to simply ready the archdat
% before the loop begins.
%
% Note that if ad is non-scalar, every element in the array will be
% readied. If the array is large, this could rapidly deplete the available
% memory.

toRelease = [ob.inmem];
[ob(toRelease).data] = deal([]);
[ob(toRelease).inmem] = deal(false);
