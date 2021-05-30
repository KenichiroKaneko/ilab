close all;
clear all;
figure()
hold on
DTHETA = pi / 3;
RR = 0.02;
capper = 1.4;
R = 2.258987252292017e-01;
Z = 3.968890906196872e-01;

Nz = 2033;
Nr = 602;
zmin = -9.985000000000001e-01;
zmax = 9.985000000000001e-01;
rmin = 1.081500000000000e-01;
rmax = 0.8880;
zz = linspace(zmin, zmax, Nz);
rr = linspace(rmin, rmax, Nr);
v = linspace(-20, 20, 21);
vars = load("../i/GSPLOT_OUTPUT/z0800_r602/merged.mat");
psi = vars.psi;
% contour(rr, zz, psi' * 1000, v);

for j = 1:6
    THETA = pi / 2.0 - DTHETA * (j - 1); % CCSでの積分は、時計回り
    RCCS(j) = R + capper * RR * cos(THETA + asin(0.0) * sin(THETA));
    ZCCS(j) = Z + RR * sin(THETA);
end

% scatter(RCCS, ZCCS, 'o')
% hold on

I = 1:3;
NCCS = 6;

RCCS(NCCS + 1) = RCCS(1);
ZCCS(NCCS + 1) = ZCCS(1);

RCCN(3 * I - 2) = (5 .* RCCS(2 * I - 1) + 5 .* RCCS(2 * I) - RCCS(2 * I + 1)) / 9;
ZCCN(3 * I - 2) = (5 .* ZCCS(2 * I - 1) + 5 .* ZCCS(2 * I) - ZCCS(2 * I + 1)) / 9;
RCCN(3 * I - 1) = RCCS(2 * I);
ZCCN(3 * I - 1) = ZCCS(2 * I);
RCCN(3 * I) = (5 .* RCCS(2 * I + 1) + 5 .* RCCS(2 * I) - RCCS(2 * I - 1)) / 9;
ZCCN(3 * I) = (5 .* ZCCS(2 * I + 1) + 5 .* ZCCS(2 * I) - ZCCS(2 * I - 1)) / 9;

RCCN = [RCCN RCCN];
ZCCN = [-ZCCN ZCCN];
scatter(RCCN, ZCCN, 'o')
xlim([0, 0.9]);
ylim([-1, 1]);
% hold on
% scatter(RCCN, -ZCCN, 'o')

%  Draw a fine curve to express the shape of CCS
MAXG = 200;
GMAX = MAXG;
DEL = 2.0 / GMAX;

count = 1;

for I = 1:3

    for J = 1:MAXG + 1
        CJM1 = J - 1;
        GII = -1.0 + DEL * CJM1;
        F1 = GII * (GII - 1.0D0) * 0.5;
        F2 = 1.0D0 - GII^2;
        F3 = GII * (GII + 1.0D0) * 0.5;
        RGI(count) = RCCS(2 * I - 1) * F1 + RCCS(2 * I) * F2 + RCCS(2 * I + 1) * F3;
        ZGI(count) = ZCCS(2 * I - 1) * F1 + ZCCS(2 * I) * F2 + ZCCS(2 * I + 1) * F3;
        count = count + 1;

    end

end

% scatter(RGI, ZGI);
