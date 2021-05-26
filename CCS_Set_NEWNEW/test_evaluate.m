% 最小２乗誤差を返すプログラム
vars = load("vars_result_600");
% vars = load("vars_result_800");
% vars = load("vars_result_1000R");
% vars = load("vars_result_2033R");

psi = vars.psi;
REF = vars.REF;
PARAM = vars.PARAM;
CCR = vars.CCR;
CCZ = vars.CCZ;

MSE = EVALUATE0(psi, REF, PARAM, CCR, CCZ)

function MSE = EVALUATE0(psi, REF, PARAM, CCR, CCZ)
    figure()
    xlim([0.1 0.9]);
    ylim([-1 1]);
    hold on

    if PARAM.CCS > 1 && abs(PARAM.Z0(1) - PARAM.Z0(2)) > 0.0
        % 2つのプラズマが離れている場合

        for i = 1:2
            j = i - 1;
            range1 = 1 + j * 50:51 + j * 50;
            range2 = 1 + j * 1016:1017 + j * 1016;

            % 再構成した磁束のLCFSを探索
            [LCFSrecon, disuseNum1, disuseNum2] = detectLCFS(psi(range1, :), CCR, CCZ(range1), -1, "reconstruction");
            % referenceのLCFSを探索
            [LCFSref, MAxisR, MAxisZ] = detectLCFS(REF.Flux(range2, :), REF.R, REF.Z(range2), 0.1656627, "reference");
            % 求めた2つの磁気面上の何点かと磁気軸からの距離をそれぞれ計算し、２乗誤差を計算
            MSE(i) = calcError(LCFSrecon, LCFSref, MAxisR, MAxisZ, CCR, CCZ(range1), REF.R, REF.Z(range2));
        end

        % [LCFSrecon, disuseNum1, disuseNum2] = detectLCFS(psi(51:101, :), CCR, CCZ(51:101), -1, "reconstruction");
        % [LCFSref, MAxisR, MAxisZ] = detectLCFS(REF.Flux(1017:2033, :), REF.R, REF.Z(1017:2033), 0.1656627, "reference");
        % MSE(2) = calcError(LCFSrecon, LCFSref, MAxisR, MAxisZ, CCR, CCZ(51:101), REF.R, REF.Z(1017:2033));
    else
        % 2つのプラズマが接近した場合

        [LCFSrecon, disuseNum1, disuseNum2] = detectLCFS(psi, CCR, CCZ, -1, "reconstruction");
        [LCFSref, MAxisR, MAxisZ] = detectLCFS(REF.Flux, REF.R, REF.Z, 0.1656627, "reference");
        MSE = calcError(LCFSrecon, LCFSref, MAxisR, MAxisZ, CCR, CCZ, REF.R, REF.Z);
    end

    MSE = sum(MSE) / length(MSE);

end

function [STRUCTURE, MAxisR, MAxisZ] = detectLCFS(psi, psiR, psiZ, limiter, type)
    psi_temp = psi;
    psi_inside = psi;

    % 磁気軸の中心を探す→refでは利用されるがreconでは要らない
    [M1 MAxisR] = min(min(psi));
    [M2 MAxisZ] = min(psi(:, MAxisR));

    if limiter == -1
        % リミターがない場合→０未満の磁束の最大値
        LCFS_min = max(max(psi_temp(psi_temp < 0))) - 0.0005;
    else
        % リミターがある場合→リミター位置の磁束
        % scatter(psiR(MAxisR), psiZ(MAxisZ), "*m");
        [m limiterR] = min(abs(psiR - limiter));
        LCFS_min = psi(MAxisZ, limiterR);
    end

    if type == "reference"
        line = "-k";
    else
        line = "-r";
    end

    psi_inside(psi_temp > LCFS_min) = 0; % LCFSの外側を０にする
    [Z, R] = ind2sub(size(psi_temp), find(psi_inside)); % LCFSの内側の座標
    j = boundary(R, Z);
    STRUCTURE.j = j;
    STRUCTURE.R = R;
    STRUCTURE.Z = Z;
    % scatter(psiR(R), psiZ(Z), 'k');

    % 他の境界プロット法
    % [B, L] = bwboundaries(psi_inside, 'noholes');

    % for k = 1:length(B)
    %     border = B{k};
    %     plot(psiR(border(:, 2)), psiZ(border(:, 1)), ':b', 'LineWidth', 1)
    % end
end

function MSE = calcError(LCFSrecon, LCFSref, MAxisR, MAxisZ, CCR, CCZ, REFR, REFZ)

    poly1 = polyshape(CCR(LCFSrecon.R(LCFSrecon.j)), CCZ(LCFSrecon.Z(LCFSrecon.j)));
    poly2 = polyshape(REFR(LCFSref.R(LCFSref.j)), REFZ(LCFSref.Z(LCFSref.j)));

    plot(poly1, 'FaceAlpha', 0.4, 'FaceColor', "c")
    plot(poly2, 'FaceAlpha', 0.4, 'FaceColor', "m")
    scatter(REFR(MAxisR), REFZ(MAxisZ));

    p = 32;
    % delta = 1.0e-8;
    delta = 0;
    l = 0.3;
    Z1 = REFZ(MAxisZ);
    R1 = REFR(MAxisR);

    for i = 1:p
        theta = 2 * (i + delta) * pi / p;

        % 磁気軸から結ぶ点
        R2 = R1 + l * sin(theta);
        Z2 = Z1 + l * cos(theta);
        % ２点を結ぶ直線
        lineseg = [R1 Z1; R2 Z2];

        % 再構成したLCFSと直線の交わり
        [in, out] = intersect(poly1, lineseg);
        x1 = in(1, 1); y1 = in(1, 2); x2 = in(end, 1); y2 = in(end, 2);
        L_recon = sqrt((x1 - x2)^2 + (y1 - y2)^2);
        % plot(in(:, 1), in(:, 2), 'b')
        % plot(out(:, 1), out(:, 2), 'r')
        % scatter(in(end, 1), in(end, 2), "ro")
        % scatter(in(1, 1), in(1, 2), "ro")

        % referenceLCFSと直線の交わり
        [in, out] = intersect(poly2, lineseg);
        x1 = in(1, 1); y1 = in(1, 2); x2 = in(end, 1); y2 = in(end, 2);
        L_ref = sqrt((x1 - x2)^2 + (y1 - y2)^2);
        % plot(in(:, 1), in(:, 2), 'b')
        % plot(out(:, 1), out(:, 2), 'r')
        % scatter(in(end, 1), in(end, 2), "yo")
        % scatter(in(1, 1), in(1, 2), "yo")

        error(i) = (L_recon - L_ref)^2;

        % plot(x, y, "k")

    end

    MSE = sum(error) / p;

end

% psiBW = psi_inside;
% psiBW(ind2sub(size(psi_inside), find(psi_inside))) = 1;
% psiBW = imfill(psiBW, 'holes');
% psiWB = abs(1 - psiBW);
% figure()
% contour(CCR, CCZ, psi, "-k", 'LevelStep', 0.0003);
% hold on
% % CCNを広げるiranaikamo
% l = 0.03;
% n = length(CCSDAT.RCCN(1, :));

% for i = 1:n
%     theta = 2 * pi / n * (i - 1) + pi / n;
%     r = l * sin(theta);
%     z = l * cos(theta);
%     CCSDAT.RCCN(:, i) = CCSDAT.RCCN(:, i) + r;
%     CCSDAT.ZCCN(:, i) = CCSDAT.ZCCN(:, i) + z;
% end
