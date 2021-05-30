close all;

A = load("GSPLOT_OUTPUT/z2033_r6022021/merged");
B = load("psiGSplotTest");
C = load("OUTPUT/z2033_r602/vars");
D = load("GSPLOT_OUTPUT/z2033_r602/merged");
psi = A.psi;
v = linspace(-20, 20, 21);
psiorig = A.psiorig0;

Nz = 2033;
Nr = 602;
zmin = -9.985000000000001e-01;
zmax = 9.985000000000001e-01;
rmin = 1.081500000000000e-01;
rmax = 0.8880;
zz = linspace(zmin, zmax, Nz);
rr = linspace(rmin, rmax, Nr);

% グラフの領域2020/12/21
vars3c = load("OUTPUT/z2033_r602/vars");
env3c = vars3c.param;
zz = linspace(env3c.zmin, env3c.zmax, env3c.Nz);
rr = linspace(env3c.rmin, env3c.rmax, env3c.Nr);

figure()
contour(rr, zz, psi' * 1000, v);

figure()
contour(rr, zz, psiorig' * 1000, v);
figure()
contour(rr, zz, B.psi' * 1000, v);

flux = A.flux;
figure()
contour(rr, zz, flux * 1000, v);

figure()
contour(rr, zz, C.psi * 1000, v);
figure()
contour(rr, zz, D.psi' * 1000, v);
