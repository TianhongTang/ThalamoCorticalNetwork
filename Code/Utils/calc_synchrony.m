% N-by-T matrix, calculate synchrony index Chi.
function chi = calc_synchrony(x)
    [N, ~] = size(x);
    if N < 1
        chi = NaN;
        return;
    end
    var_ave = var(mean(x, 1), 0, 2);
    var_ind = mean(var(x, 0, 2), 1);
    chi = sqrt(var_ave / var_ind);
end