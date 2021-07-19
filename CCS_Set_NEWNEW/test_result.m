close all
v = load("vars_result_ima1ji");
REF = v.REF;
CCR = v.CCR;
CCZ = v.CCZ;
psi = v.psi;

Flux = REF.Flux - 0.1;
% psi2 = psi + 0.00002;

figure
contour(CCR, CCZ, psi2, '-k', 'LevelStep', 0.0003);
hold on
contour(REF.R, REF.Z, REF.Flux, '--m', 'LevelStep', 0.0003); % 正解

contour(CCR, CCZ, psi2, [0 0], 'LineColor', 'c', 'LineWidth', 2);
contour(REF.R, REF.Z, REF.Flux, [0 0], 'LineColor', 'm', 'LineWidth', 2);

hold off
xlabel({'r (m)'});
ylabel({'z (m)'});
title("Reconstructed flux")
axis equal

%% Data positions for 2D plot
PLOT.R = 0.0:0.001:1;
PLOT.Z = -1:0.002:1;

figure
v = linspace(-20, 20, 41);
contour(CCR, CCZ, psi * 1000, v, '-k');
hold on
Flux = REF.Flux;
% Flux = Flux(:, 1:602);
Flux = imresize(Flux, [length(CCZ), length(CCR)]);
contour(CCR, CCZ, Flux * 1000, v, '--m'); % ????
hold off
xlabel({'r (m)'});
ylabel({'z (m)'});
title("Reconstructed flux")
axis equal

%% Data positions for 2D plot
PLOT.R = 0.0:0.001:1;
PLOT.Z = -1:0.002:1;
