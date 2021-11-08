vars = load('vars_result_actualsensorpoosition');
figure
hold on
contour(vars.CCR, vars.CCZ, vars.psi * 1000, vars.v, '-k');
contour(vars.REF.R, vars.REF.Z, vars.REF.Flux * 1000, vars.v, '--m'); % 正解
contour(vars.CCR, vars.CCZ, vars.psi, [0 0], 'LineColor', 'c', 'LineWidth', 2);
contour(vars.REF.R, vars.REF.Z, vars.REF.Flux, [0 0], 'LineColor', 'm', 'LineWidth', 2);
plot(vars.VV(:, 1), vars.VV(:, 2), '-k'); % 容器壁
legend('Reconstructed Flux', 'Reference Flux', 'Reconstructed LCFS', 'Reference LCFS', 'Vacuum vessel wall');
hold off
xlabel({'r (m)'});
ylabel({'z (m)'});
title("Reconstructed flux")
axis equal
