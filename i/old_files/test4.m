psi602 = load('psi602', 'psi').psi;
psi800 = load('psi800', 'psi').psi;

size(psi800)

padding = zeros(2033, 198);

psi602 = [psi602 padding];

z = linspace(param.zmin, param.zmax, param.Nz);
r = linspace(param.rmin, param.rmax, param.Nr);
v = linspace(-50, 50, 21);

figure()
contour(r, z, (psi800-psi602) * 1000, v)


A = psi800 - psi602;

save("psiOutside.mat", "A");