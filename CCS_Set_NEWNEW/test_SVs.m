close all;

SVs = load('test_SVs').SVS;

figure()
plot(1:87, log(SV));

figure()
plot(log(SV))

f = log(SV);
x = 1:87;
y = -x;

KUP = curvature(SVs)

