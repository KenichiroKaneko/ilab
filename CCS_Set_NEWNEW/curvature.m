function KUP = CURVATURE(X, Y);

    x = log10(X);
    y = log10(Y);

    % 前後の無視する数/点群の数、num < cut
    cut = 5;
    num = 1;

    len = length(X);
    m = 0;
    KUP = 1;

    % ３点で求める曲率
    for i = 1:len - cut * 2

        j = i + cut;
        p0 = [x(j - 1); y(j - 1)];
        p1 = [x(j); y(j)];
        p2 = [x(j + 1); y(j + 1)];

        A = p2 - p1;
        B = p0 - p1;

        Numerator = 2 * abs(det([p2 - p1 p0 - p1]));
        Denominator = norm(p2 - p1) * norm(p0 - p1) * norm(p2 - p0);

        K(i) = Numerator / Denominator;

    end

    % 傾きの差の最大値
    for i = 1:length(x) - cut * 2

        j = i + cut;
        Left = (y(j) - y(j - 1)) / (x(j) - x(j - 1));
        Right = (y(j + 1) - y(j)) / (x(j + 1) - x(j));

        tilt(i) = Right - Left;

        if m < tilt(i)
            m = tilt(i);
            KUP = i;
        end

    end

    % num点で求める曲率
    for i = 1:length(x) - cut * 2

        j = i + cut;
        xs = x(j - num:j + num);
        ys = y(j - num:j + num);

        [cx, cy, r] = CircleFitting(xs, ys);

        R(i) = r;

    end

    % y=x-1000と一番近い点 ax + by + c = 0
    y_intercept = 10;
    a = 1;
    b = 1;
    c = y_intercept;

    for i = 1:length(x) - cut * 2
        j = i + cut;
        d(i) = abs(a * x(j) + b * y(j) + c) / norm([a b]);
    end

    [Mk, Ik] = max(K)
    Ik = Ik + cut;
    KUPIndex_K = len + 1 - Ik;

    KUP = KUP + cut;

    [Mr, Ir] = max(R)
    I = Ir + cut;
    KUPIndex = len + 1 - I;

    [My Iy] = min(d);
    Iy = Iy + cut;
    KUPIndex_Y = len + 1 - Iy;
    abcx = -10:0;
    abcy = -(abcx) - y_intercept;

    pad = zeros(1, cut);
    K = [pad K pad];
    tilt = [pad tilt pad];
    R = [pad R pad];
    d = [pad d pad];

    figure()
    % ｘーｙと曲率最大のポイント赤色
    loglog(X, Y, 's')
    hold on
    scatter(X(KUPIndex_K), Y(KUPIndex_K), "r*")
    scatter(X(KUP), Y(KUP), 'b*')
    scatter(X(KUPIndex), Y(KUPIndex), "g*")
    scatter(X(Iy), Y(Iy), "y*")
    hold off
    % 曲率、傾きをｘーｙに合わせてプロットしたもの
    figure()
    hold on
    % plot(x, y * 10, "black")
    plot(x, K, "r")
    plot(x, tilt, "b")
    plot(x, R, "g")
    plot(x, d, "y");
    scatter(x(Ik), K(Ik), "r*")
    scatter(x(KUP), tilt(KUP), "b*")
    scatter(x(I), R(I), 'g*')
    scatter(x(Iy), R(Iy), 'y*')
    hold off
    % 曲率、傾きをプロットしたもの
    figure()
    hold on
    % x = fliplr(x);
    % plot(y * 10, "black")
    plot(K, "r")
    plot(tilt, "b")
    plot(R, "g")
    plot(d, "y")
    scatter(Ik, K(Ik), "r*");
    scatter(KUP, tilt(KUP), "b*");
    scatter(I, R(I), "g*");
    scatter(Iy, R(Iy), "y*");
    hold off

    figure()
    scatter(x, y);
    hold on
    plot(abcx, abcy)
    scatter(x(Iy), y(Iy), "r*")
    scatter(x(KUPIndex_Y), y(KUPIndex_Y), "b*")

    figure()
    hold on
    scatter(1:length(d), d);
    scatter(Iy, d(Iy), "r*")
    scatter(KUPIndex_Y, d(KUPIndex_Y), "b*")

    KUP = Iy;

end

function [cx, cy, r] = CircleFitting(x, y)
    %CIRCLEFITTING 最小二乗法による円フィッテングをする関数
    % input: x,y 円フィッティングする点群
    % output cx 中心x座標
    %        cy 中心y座標
    %        r  半径
    % 参考
    % 一般式による最小二乗法（円の最小二乗法）　画像処理ソリューション
    % http://imagingsolution.blog107.fc2.com/blog-entry-16.html

    sumx = sum(x);
    sumy = sum(y);
    sumx2 = sum(x.^2);
    sumy2 = sum(y.^2);
    sumxy = sum(x .* y);

    F = [sumx2 sumxy sumx;
        sumxy sumy2 sumy;
        sumx sumy length(x)];

    G = [-sum(x.^3 + x .* y.^2);
        -sum(x.^2 .* y + y.^3);
        -sum(x.^2 + y.^2)];

    T = F \ G;

    cx = T(1) / -2;
    cy = T(2) / -2;
    r = sqrt(cx^2 + cy^2 - T(3));

end
