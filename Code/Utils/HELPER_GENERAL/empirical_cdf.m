function [y] = empirical_cdf(stat,ref)
% [y] = empirical_cdf(stat,ref)
% cdf of the empirical distribution specified by 'ref'
%
% ways of calling:
%
% ref as a data vector, used to compute the cdf.
% this approximates normcdf(x,0,1):
% > ref = randn(1,1000);
% > y = empirical_cdf(x,ref);
%
% ref as a data vector for each entry in 
% this approximates poisscdf(x,1:50):
% > ref = poissrnd(repmat(1:50,1000,1));
% > y = empirical_cdf(x,ref)
%

if isnumeric(ref)
    if isvector(ref)
        ncols = numel(stat);
        stat_collapse = reshape(stat,1,numel(stat));
        y = zeros(size(stat_collapse));
        for i = 1:length(stat_collapse)
            y(i) = mean(ref <= stat_collapse(i));
        end;
        y = reshape(y,size(stat));
    else
        assert(isvector(stat) && ndims(ref) == 2);
        if size(stat,1) == 1
            assert(size(stat,2) == size(ref,2));
            stat_big = repmat(stat,size(ref,1),1);
            y = mean(ref <= stat_big,1);
        else
            assert(size(stat,1) == size(ref,1));
            stat_big = repmat(stat,1,size(ref,2));
            y = mean(ref <= stat_big,2);
        end;
    end;
elseif iscell(ref)
    stat_collapse = reshape(stat,1,numel(stat));
    count = zeros(size(stat_collapse));
    for r = 1:length(ref)
        ref_collapse = reshape(ref{r},size(stat_collapse));
        count = count + (ref_collapse <= stat_collapse);
    end;
    y = count ./ length(ref);
    y = reshape(y,size(stat));
else
    error('unknown input type');
end;
