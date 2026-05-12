function [heart_rate, TI_flg] = ECGSigProcFB(signal)

TI_flg = false;
samp_rate = size(signal,1) / (signal(end,1)-signal(1,1));
clf('reset'); subplot(2, 1, 1); plot(signal(:,1),signal(:,2));

QRS=nqrsdetect(signal(:,2), samp_rate);

hold on; plot(signal(QRS,1),signal(QRS,2),'r*');
title('ECG Signal');

m = size(QRS, 1);

if m >= 2
    rr_sig = signal(QRS(1):QRS(2), :);
    subplot(2, 1, 2); plot(rr_sig(:, 1), rr_sig(:, 2));
    title('R-R Interval Signal');


    % T-wave inversion detection.
    min_idx = find(rr_sig(:, 2) == min(rr_sig(:, 2)));

    % Position is not within 10% of R-wave.
    if (rr_sig(min_idx, 1) - rr_sig(1, 1)) >= (0.1 * (rr_sig(end, 1) - rr_sig(1, 1))) && (rr_sig(end, 1) - rr_sig(min_idx, 1)) >= (0.1 * (rr_sig(end, 1) - rr_sig(1, 1)))
        hold on; plot(rr_sig(min_idx, 1), rr_sig(min_idx, 2), 'r*');
        TI_flg = true;
    end
end

if m >= 2
    heart_rate = 60 / mean(diff(signal(QRS,1)));
else
    heart_rate = 0;
end

hold off;

