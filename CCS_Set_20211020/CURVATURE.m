function KUP = CURVATURE(X, Y);

    x = log10(X);
    y = log10(Y);

    % 前後の無視する数/点群の数、num < cut
    cut = 10;
    num = 1;

    len = length(X);
    m = 0;
    KUP = 1;

    % y=x-10と一番近い点 ax + by + c = 0
    y_intercept = 10;
    a = 1;
    b = 1;
    c = y_intercept;

    for i = 1:length(x) - cut * 2
        j = i + cut;
        d(i) = abs(a * x(j) + b * y(j) + c) / norm([a b]);
    end

    [My Iy] = min(d);
    Iy = Iy + cut;
    KUPIndex_Y = len + 1 - Iy;

    pad = zeros(1, cut);
    d = [pad d pad];

    % ｘーｙと曲率最大のポイント赤色
    % figure()
    % loglog(X, Y, 's')
    % hold on
    % scatter(X(Iy), Y(Iy), "y*")
    % hold off
    KUP = Iy;

end
