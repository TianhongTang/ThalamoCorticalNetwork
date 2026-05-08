function plot_csd(csd_in, time_vec, tlim, ax)

csd_in = [csd_in(1,:) ; csd_in ; csd_in(end,:)];
nele = size(csd_in,1);
csd_in = H_2DSMOOTH(csd_in);

limi = nanmax(nanmax(abs(csd_in)));

imagesc(time_vec, 1:size(csd_in,1), csd_in);

set(gca,'xlim', tlim, 'ydir', 'rev', ...
    'ylim', [min(find(~isnan(csd_in(:,1)))) max(find(~isnan(csd_in(:,1))))])

end