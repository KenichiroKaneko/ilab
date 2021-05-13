function KUP = curvature(SVs)
    SV = SVs(:, 1);
    f = log(SV);

    max = 0;
    KUP = 0;

    for i = 10:(length(f) - 10)
        Left = f(i) - f(i - 1);
        Right = f(i + 1) - f(i);

        if max < (Right - Left)
            max = (Right - Left);
            KUP = i;
        end

    end

end
