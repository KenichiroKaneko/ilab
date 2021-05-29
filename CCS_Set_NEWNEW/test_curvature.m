close all;
clear all;
format longG
vars = load("vars_inLcurve");

X = vars.residual;
Y = vars.P_norm;
x = log10(X);
y = log10(Y);
% X = 0:0.1:10;
% Y = sin(X) + X.^2;
% x = X;
% y = Y;
f = fit(x', y', "cubicspline");
fnplt(f.p, "-k", 1)
hold on
scatter(x, y);
hold off
% 前後の無視する数/点群の数、num < cut
cut = 5;
num = 1;
pad = zeros(1, cut);

len = length(X);
m = 0;
KbyTilt = 1;

% １
% ３点で求める曲率
if 1
    spX = x(cut):(x(end - cut) - x(cut)) / 100:x(end - cut);
    spY = spline(x, y, spX);

    for i = 1:100 - cut * 2
        j = i + num;
        xs = spX(j - num:j + num);
        ys = spY(j - num:j + num);
        [cx, cy, r] = CircleFitting(xs, ys);
        spR(i) = 1 / r;
    end

    figure()
    plot(spX, spY)
    hold on
    % figure()
    % plot(1:100 - cut * 2, spR)

    x = spX;
    y = spY;

    % for i = 1:len - cut * 2
    for i = 1:100 - cut * 2
        j = i + cut;
        p0 = [x(j - 1); y(j - 1)];
        p1 = [x(j); y(j)];
        p2 = [x(j + 1); y(j + 1)];

        A = p2 - p1;
        B = p0 - p1;

        Numerator = 2 * abs(det([p2 - p1 p0 - p1]));
        Denominator = norm(p2 - p1) * norm(p0 - p1) * norm(p2 - p0);

        Kby3(i) = Numerator / Denominator;

    end

    [Mk, IndKby3] = max(Kby3)
    IndKby3 = IndKby3 + cut;
    KUPIndex_K = len + 1 - IndKby3;
    KUPIndex_K = IndKby3;

    scatter(spX(KUPIndex_K), spY(KUPIndex_K))

    % figure('Name', '３点で求める曲率')
    % loglog(X, Y, 's')
    % hold on
    % scatter(X(KUPIndex_K), Y(KUPIndex_K), "r*")
    % Kby3 = [pad Kby3 pad];
end

% ２
% 傾きの差の最大値
if 0

    for i = 1:length(x) - cut * 2

        j = i + cut;
        Left = (y(j) - y(j - 1)) / (x(j) - x(j - 1));
        Right = (y(j + 1) - y(j)) / (x(j + 1) - x(j));

        tilt(i) = Right - Left;

        if m < tilt(i)
            m = tilt(i);
            KbyTilt = i;
        end

    end

    KbyTilt = KbyTilt + cut;
    figure('Name', '前後の傾きの最大値法')
    loglog(X, Y, 's')
    hold on
    scatter(X(KbyTilt), Y(KbyTilt), 'b*')
    % tilt = [pad tilt pad];
    % scatter(x(KbyTilt), tilt(KbyTilt), "b*")
end

% ３
% num点で求める曲率
if 0

    for i = 1:length(x) - cut * 2

        j = i + cut;
        xs = x(j - num:j + num);
        ys = y(j - num:j + num);

        [cx, cy, r] = CircleFitting(xs, ys);

        R(i) = 1 / r;

    end

    [Mr, Ir] = max(R)
    I = Ir + cut;
    KUPIndex = len + 1 - I;
    figure('Name', 'ネットで見つけた曲率求めるやつ①')
    loglog(X, Y, 's')
    hold on
    scatter(X(KUPIndex), Y(KUPIndex), "g*")
    R = [pad R pad];
end

% ４
% 曲線をスプライン近似してやる
if 1
    spX = x(cut):(x(end - cut) - x(cut)) / 100:x(end - cut);
    spY = spline(x, y, spX);

    for i = 1:100 - cut * 2
        j = i + num;
        xs = spX(j - num:j + num);
        ys = spY(j - num:j + num);
        [cx, cy, r] = CircleFitting(xs, ys);
        spR(i) = 1 / r;
    end

    figure()
    plot(spX, spY)
    figure()
    plot(1:100 - cut * 2, spR)
end

% ５
% y=x-1000と一番近い点 ax + by + c = 0
if 0
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
    abcx = -10:0;
    abcy = -(abcx) - y_intercept;
    d = [pad d pad];

    figure('Name', '直線との距離y=-x-10')
    loglog(X, Y, 's')
    hold on
    scatter(X(Iy), Y(Iy), "y*")
    hold off
end

% 曲率、傾きをｘーｙに合わせてプロットしたもの
% figure()
% hold on
% % plot(x, y * 10, "black")
% plot(x, Kby3, "r")
% plot(x, tilt, "b")
% plot(x, R, "g")
% plot(x, d, "y");
% scatter(x(IndKby3), Kby3(IndKby3), "r*")
% scatter(x(I), R(I), 'g*')
% scatter(x(Iy), R(Iy), 'y*')
% hold off
% % 曲率、傾きをプロットしたもの
% figure()
% hold on
% % x = fliplr(x);
% % plot(y * 10, "black")
% plot(Kby3, "r")
% plot(tilt, "b")
% plot(R, "g")
% plot(d, "y")
% scatter(IndKby3, Kby3(IndKby3), "r*");
% scatter(KbyTilt, tilt(KbyTilt), "b*");
% scatter(I, R(I), "g*");
% scatter(Iy, R(Iy), "y*");
% hold off

% function KUP = curvature(f)

%     f;

%     max = 0;
%     KUP = 0;

%     for i = 5:(length(f) - 5)
%         Left = f(i) - f(i - 1);
%         Right = f(i + 1) - f(i);

%         if max < (Right - Left)
%             max = (Right - Left);
%             KUP = i;
%         end

%     end

% end

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
