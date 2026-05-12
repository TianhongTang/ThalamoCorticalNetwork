function [y] = shuffle(x,mode,n,groups,nanmode)
% [y] = shuffle(x,mode,[n,groups,nanmode])
%  shuffle a matrix. basically a souped-up version of randperm.
%
% inputs:
%  x - the matrix to be shuffled
%
%  mode - there are five modes:
%   'rows' - shuffle the rows
%   'cols' - shuffle the cols
%   'inrows' - shuffle the elements within each row
%   'incols' - shuffle the elements within each column
%   'all' - shuffle all elements
%  (default: 'all')
%
%  n - the number of shufflings.
%       if x is a (1xN) vector, returns the shufflings as a (nxN) matrix
%       if x is a (Nx1) vector, returns the shufflings as a (Nxn) matrix
%       if x is a matrix, returns them as a (nx1) cell array of shuffled matrices
%      you can also specify n as the string 'one'. 
%       If so, then it will do a single shuffling and will return it 
%       without wrapping it in a cell array, even if the result is a matrix. 
%      (default = 'one')
%
%  groups - matrix assigning a 'group number' to each shuffled item.
%   Only items with the same group number may be interchanged during the
%   shuffling process. For each mode, groups must have a specific size:
%       'rows' - [size(x,1) 1]
%       'cols' - [1 size(x,2)]
%       'inrows' - size(x) (or [1 size(x,2)] will be replicated to fit size(x))
%       'incols' - size(x) (or [size(x,1) 1] will be replicated to fit size(x))
%       'all' - size(x)
%   (default = a matrix of ones. That is, all items are put in a single group 
%    thus allowing all possible shufflings)
%
%  nanmode - how to handle NaN entries.
%   'ignore' - do not change the position of NaN entries
%   'shuffle' - shuffle them just like any other entries
%   (default: 'shuffle')
%
% NOTE: changed 2010-02-10 to replace calls to 'randperm' with
%  hardcoded version of the function ([jnk,rperm] = sort(rand(1,groupsize));)
%  which greatly reduces running time b/c don't need so many function
%  calls.
%
% examples:
%
% > test = [1 2 3; ...
%           4 5 6; ...
%           7 8 9];
%
% to shuffle all elements:
% > shuffle(test)
%     ans =
% 
%          2     3     8
%          6     1     4
%          5     7     9
%
% to shuffle the rows:
% > shuffle(test,'rows')
%     ans =
%          1     2     3
%          7     8     9
%          4     5     6
%
% to shuffle separately within each row:
% > shuffle(test,'inrows')
%     ans =
%          2     3     1
%          5     4     6
%          9     7     8
%
% to shuffle separately the even and odd numbers:
% > shuffle(test,'all','one',mod(test,2))
%     ans =
%          9     6     1
%          4     3     8
%          5     2     7
%
% to shuffle separately within each column,
%  only interchanging pairs of elements which have
%  the same oddness
% > shuffle(test,'incols','one',mod(test,2))
%     ans =
%          7     8     3
%          4     5     6
%          1     2     9

if nargin < 2 || isempty(mode)
    mode = 'all';
end;
if nargin < 3 || isempty(n)
    n = 'one';
end;
if nargin < 4 || isempty(groups)
    switch mode
        case {'rows','incols'}
            groups = ones(size(x,1),1);
        case {'cols','inrows'}
            groups = ones(1,size(x,2));
        case 'all'
            groups = ones(size(x));
        otherwise
            error('unknown mode');
    end;
end;
if nargin < 5 || isempty(nanmode)
    nanmode = 'shuffle';
end;

do_not_wrap_in_cell_array = false;
if ischar(n)
    do_not_wrap_in_cell_array = true;
    n = 1;
end;

if strcmp(mode,'all')
    assert(all(size(groups) == size(x)),'bad group size');
else
    assert(ndims(x) == 2,['mode ' mode ' implemented only for 2d matrices (try ''all'')']);
    if strcmp(mode,'rows')
        assert(all(size(groups) == [size(x,1) 1]),'bad group size');
    elseif strcmp(mode,'cols')
        assert(all(size(groups) == [1 size(x,2)]),'bad group size');
    elseif strcmp(mode,'inrows')
        if all(size(groups) == size(x))
        elseif all(size(groups) == [1 size(x,2)])
            groups = repmat(groups,size(x,1),1);
        else
            error('bad group size');
        end;
    elseif strcmp(mode,'incols')
        if all(size(groups) == size(x))
        elseif all(size(groups) == [size(x,1) 1])
            groups = repmat(groups,1,size(x,2));
        else
            error('bad group size');
        end;
    else
        error('unknown mode');
    end;
end

assert(~any(isnan(groups(:))),'group ids may not be NaN!');

if strcmp(mode,'all')
    % collapse x to a row vector, shuffle it, then un-collapse it
    x_oldsize = size(x);
    was_rowvector = isrowvector(x);
    was_colvector = iscolvector(x);
    
    x = reshape(x,1,numel(x));
    groups = reshape(groups,1,numel(groups));
    ytemp = shuffle(x,'cols',n,groups,nanmode);
    
    if was_rowvector
        y = ytemp;
    elseif was_colvector
        y = ytemp';
    else
        y = cell(n,1);
        for s = 1:n
            y{s} = reshape(ytemp(s,:),x_oldsize);
        end;
    end;
else
    unique_groups = unique(groups);
    
    % do we ignore NaNs? If so, put them in a special 'NaN' group
    %  which will be ignored by all our other methods
    if strcmp(nanmode,'ignore')
        xnans = isnan(x);
        if any(xnans(:))
            if islogical(groups)
                groups = double(groups);
            end;
            groups(xnans) = nan;
        end;
        % initialize output to NaN
        if isrowvector(x)
            y = nans(n,size(x,2));
        elseif iscolvector(x)
            y = nans(size(x,1),n);
        else
            y = cell(n,1);
            for s = 1:n
                y{s} = nans(size(x));
            end;
        end;
    else
        % for speed, initialize output to zeros, instead of NaNs
        if isrowvector(x)
            y = zeros(n,size(x,2));
        elseif iscolvector(x)
            y = zeros(size(x,1),n);
        else
            y = cell(n,1);
            for s = 1:n
                y{s} = zeros(size(x));
            end;
        end;
    end;
    
    for u = 1:length(unique_groups)
        curgroup = groups == unique_groups(u);
        
        if strcmp(mode,'rows') || strcmp(mode,'incols')
            if isrowvector(x)
                % shuffling within columns has no effect on a row vector
                for s = 1:n
                    y(s,curgroup) = x(curgroup);
                end;
            elseif iscolvector(x)
                % repeatedly shuffle the single column
                curgroup_vals = x(curgroup);
                groupsize = length(curgroup_vals);
                for s = 1:n
                    [jnk,rperm] = sort(rand(1,groupsize));
                    y(curgroup,s) = curgroup_vals(rperm);
                end;
            else
                if strcmp(mode,'incols')
                    % for each column, find elements in the current group,
                    %  and repeatedly shuffle their values
                    for c = 1:size(x,2)
                        curgroup_col = curgroup(:,c);
                        curgroup_vals = x(curgroup_col,c);
                        groupsize = size(curgroup_vals,1);
                        for s = 1:n
                            [jnk,rperm] = sort(rand(1,groupsize));
                            y{s}(curgroup_col,c) = curgroup_vals(rperm);
                        end;
                    end;
                else
                    % repeatedly shuffle the rows
                    curgroup_vals = x(curgroup,:);
                    groupsize = size(curgroup_vals,1);
                    for s = 1:n
                        [jnk,rperm] = sort(rand(1,groupsize));
                        y{s}(curgroup,:) = curgroup_vals(rperm,:);
                    end;
                end;
            end;
        elseif strcmp(mode,'cols') || strcmp(mode,'inrows')
            if isrowvector(x)
                % repeatedly shuffle the single row
                curgroup_vals = x(curgroup);
                groupsize = length(curgroup_vals);
                for s = 1:n
                    [jnk,rperm] = sort(rand(1,groupsize));
                    y(s,curgroup) = curgroup_vals(rperm);
                end;
            elseif iscolvector(x)
                % shuffling within rows has no effect on a column vector
                for s = 1:n
                    y(curgroup,s) = x(curgroup);
                end;
            else
                if strcmp(mode,'inrows')
                    % for each row, find elements in the current group,
                    %  and repeatedly shuffle their values
                    for r = 1:size(x,1)
                        curgroup_row = curgroup(r,:);
                        curgroup_vals = x(r,curgroup_row);
                        groupsize = size(curgroup_vals,2);
                        for s = 1:n
                            [jnk,rperm] = sort(rand(1,groupsize));
                            y{s}(r,curgroup_row) = curgroup_vals(rperm);
                        end;
                    end;
                else
                    % repeatedly shuffle the columns
                    curgroup_vals = x(:,curgroup);
                    groupsize = size(curgroup_vals,2);
                    for s = 1:n
                        [jnk,rperm] = sort(rand(1,groupsize));
                        y{s}(:,curgroup) = curgroup_vals(:,rperm);
                    end;
                end;
            end;
        else
            error('unknown mode');
        end;
    end;
end;

if do_not_wrap_in_cell_array && iscell(y) && (length(y) == 1)
    y = y{1};
end;