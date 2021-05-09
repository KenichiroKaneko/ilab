figure()
DTHETA = pi / 3;
RR = 0.03;
capper = 1.4;
R = 0.22;
Z = 0.38;

for j = 1:6

    THETA = pi / 2.0 - DTHETA * (j - 1); % CCSでの積分は、時計回り
    RCCS(j) = R + RR * cos(THETA + asin(0) * sin(THETA));
    ZCCS(j) = Z + capper * RR * sin(THETA);
end

scatter(RCCS, ZCCS, 'o')
hold on

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


% polarscatter(cos(THETA + 0 * sin(THETA)), j)
scatter(RCCN, ZCCN, '*')
% xlim([0, 0.9]);
% ylim([-1, 1]);
hold on

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

scatter(RGI, ZGI);