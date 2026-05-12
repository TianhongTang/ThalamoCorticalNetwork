function [o_out, p_out] = D_CSDBAND_BASIC(data_in, fs, spc, band)

% data_in: 2D (channel x sample) or 3D (channel x sample x trial) matrix of data in uV
% fs: sampling frequency in Hz (e.g., 1000)
% spc: spacing of electrode contacts in mm (e.g., 0.1)
% band in Hz (e.g., [15 30])

% o_out: oscilatory current activity in nA/mm3
% p_out: power of current oscillation in nA/mm3

cndt = 0.0004;

lpc1 = band(2);
hpc1 = band(1);
lpc2 = band(1)/2;

csd_t = data_in;
nChan = size( csd_t, 1 ) * spc;
dChan = spc : spc : nChan;

nE = length( dChan );
d = mean( diff( dChan ) );

t_csd = [];
for i = 1 : nE - 2
    for j = 1 : nE
        if i == (j - 1); t_csd( i, j ) = -2 / d^2;
        elseif abs( i - j + 1) == 1; t_csd( i, j ) = 1 / d^2;
        else; t_csd( i, j ) = 0;
        end
    end
end

o_out = nan(size(csd_t,1)-2, size(csd_t,2), size(csd_t,3));

for i = 1 : size(data_in, 3)
    o_out(:,:,i) = cndt .* (t_csd * csd_t(:,:,i)) .* 1000;
end

if size(o_out, 3) > 1
    was3D = true;
    dims = [size(o_out,1) size(o_out,2) size(o_out,3)];
    o_out = permute(o_out, [2 1 3]);
    o_out = reshape(o_out, [size(o_out,1), ...
        size(o_out,2)*size(o_out,3)]);
else
    was3D = false;
    o_out = o_out';
end

hWn = hpc1 / (fs/2);
[ bwb, bwa ] = butter( 4, hWn, 'high' );
o_out = filtfilt( bwb, bwa, o_out);

lWn = lpc1 / (fs/2);
[ bwb, bwa ] = butter( 4, lWn, 'low' );
o_out = filtfilt( bwb, bwa, o_out );

p_out = abs(o_out);

lWn = lpc2 / (fs/2);
[ bwb, bwa ] = butter( 4, lWn, 'low' );
p_out = filtfilt( bwb, bwa, p_out );

if was3D
    o_out = reshape(o_out, [dims(2), dims(1), dims(3)]);
    o_out = permute(o_out, [2 1 3]);
    p_out = reshape(p_out, [dims(2), dims(1), dims(3)]);
    p_out = permute(p_out, [2 1 3]);
else
    o_out = o_out';
    p_out = p_out';
end

o_out = cat(1,nan(1,size(o_out,2),size(o_out,3)), ...
    o_out, ...
    nan(1,size(o_out,2),size(o_out,3)));

p_out = cat(1,nan(1,size(p_out,2),size(p_out,3)), ...
    p_out, ...
    nan(1,size(p_out,2),size(p_out,3)));

end