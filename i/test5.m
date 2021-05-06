x = [1, 3, 5, 6, 5, 2, 1];
y = [1, 2, 3, 4, 5, 6, 7];
[x, y] = prepareCurveData(x, y);
f = fit(x, y, 'smoothingspline');
figure()
plot(f, x, y)
