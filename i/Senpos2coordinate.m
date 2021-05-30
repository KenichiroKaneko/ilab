% UTST内のセンサー位置を203x602のメッシュ座標に変換する

f = dlmread("SENPOS0.txt");
o = 'sensorCoordinate0.txt';
r = f(:, 1);
z = f(:, 2);

zmin = -9.985000000000001e-01;
zmax = 9.985000000000001e-01;
rmin = 1.081500000000000e-01;
rmax = 6.940000000000000e-01;
delr = 9.747920133111480e-04;
delz = 9.827755905511811e-04;

r = r - rmin;
z = z - zmin;

r = round(r/delr) + 1;
z = round(z/delz);

A = [r z];

writematrix(A, o);
