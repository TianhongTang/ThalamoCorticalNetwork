function histmat  = hist2(r, c, redges, cedges, z)
% function histmat  = hist2(r, c, redges, cedges, z)
%
% Extract 2D histogram data containing the number of events
% of [x , y] pairs that fall in each bin of the grid defined by 
% xedges and yedges. The edges are vectors with monotonically 
% non-decreasing values.  
%
%EXAMPLE 
%
% events = 1000000;
% x1 = sqrt(0.05)*randn(events,1)-0.5; x2 = sqrt(0.05)*randn(events,1)+0.5;
% y1 = sqrt(0.05)*randn(events,1)+0.5; y2 = sqrt(0.05)*randn(events,1)-0.5;
% x= [x1;x2]; y = [y1;y2];
%
%For linearly spaced edges:
% xedges = linspace(-1,1,64); yedges = linspace(-1,1,64);
% histmat = hist2(x, y, xedges, yedges);
% figure; pcolor(xedges,yedges,histmat'); colorbar ; axis square tight ;
%
%For nonlinearly spaced edges:
% xedges_ = logspace(0,log10(3),64)-2; yedges_ = linspace(-1,1,64);
% histmat_ = hist2(x, y, xedges_, yedges_);
% figure; pcolor(xedges_,yedges_,histmat_'); colorbar ; axis square tight ;

% University of Debrecen, PET Center/Laszlo Balkay 2006
% email: balkay@pet.dote.hu
%
%NOTE from Ethan: I stopped the output matrix from being transposed.
% before: histmat(1,3) ~ first y-bin and the third x-bin.
%  after: histmat(1,3) ~ first x-bin and the third y-bin.
% 
% so the result histmat is a numel(xedges) x numel(yedges) matrix.

do_mean_instead_of_sum = false;

if nargin < 4
    error ('At least four input arguments are required!');
    return;
end
assert(isequal(size(r),size(c)),'The size of r and c should be same!');
if nargin >= 5
    assert(isequal(size(r),size(z)),'The size of r and c and z should be same!');
    do_mean_instead_of_sum = true;
end;


r = r(:);
c = c(:);

[rn, rbin] = histc(r,redges);
[cn, cbin] = histc(c,cedges);

out_of_range = (rbin==0) | (cbin==0);
rbin = rbin(~out_of_range);
cbin = cbin(~out_of_range);

if do_mean_instead_of_sum
    z = z(:);
    z = z(~out_of_range);
    hist_n = zeros(numel(redges),numel(cedges));
    hist_sum = zeros(numel(redges),numel(cedges));
    for i = 1:numel(rbin)
        hist_n(rbin(i),cbin(i)) = hist_n(rbin(i),cbin(i))+1;
        hist_sum(rbin(i),cbin(i)) = hist_sum(rbin(i),cbin(i)) + z(i);
    end;

    histmat = hist_sum ./ hist_n;
else
    histmat = zeros(numel(redges),numel(cedges));
    for i = 1:numel(rbin)
        histmat(rbin(i),cbin(i)) = histmat(rbin(i),cbin(i))+1;
    end;
end;


% 
% % OLD CODE
% % %xbin, ybin zero for out of range values 
% % % (see the help of histc) force this event to the 
% % % first bins
% % xbin(find(xbin == 0)) = 1;
% % ybin(find(ybin == 0)) = 1;
% 
% % ignore out-of-range values
% out_of_range = (xbin==0) | (ybin==0);
% xbin = xbin(~out_of_range);
% ybin = ybin(~out_of_range);
% 
% xnbin = length(xedges);
% ynbin = length(yedges);
% 
% if xnbin >= ynbin
%     xy = ybin*(xnbin) + xbin;
%       indexshift =  xnbin; 
% else
%     xy = xbin*(ynbin) + ybin;
%       indexshift =  ynbin; 
% end
% 
% %[xyuni, m, n] = unique(xy);
% xyuni = unique(xy);
% hstres = histc(xy,xyuni);
% histmat = zeros(xnbin,ynbin);
% histmat(xyuni-indexshift) = hstres;
% 
% % commented out by ethan
% % histmat = histmat';
