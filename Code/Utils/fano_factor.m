function FF = fano_factor(X)
    % X: (neurons, time)
    [N, ~] = size(X);
    FF = zeros(N, 1);

    for i = 1:N
        isi = diff(find(X(i, :)>0));
        if length(isi)>1
            FF(i) = var(isi) / (mean(isi)^2);
        else
            FF(i) = NaN;
        end
    end
end