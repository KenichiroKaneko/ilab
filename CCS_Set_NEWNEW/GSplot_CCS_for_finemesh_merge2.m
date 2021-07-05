function GSplot_CCS_for_finemesh_merge2(dirname)

    % CCS03cはUTST装置全体の解
    % CCS05は（途中に仮想壁を置いた）部分的な解
    % 部分的な解を全体に拡張する際に、CCS03cのメッシュデータが必要
    % fluxは平衡解、flux0はpsiの初期状態
    close all;

    saveflag = "modan"; % origin or modan or 0
    dispFigure = 0;
    realFlag = 0;
    % 保存するdir
    id = extractBetween(dirname, 2, 5);
    save_dir = "CCS_dataset_gs/UTST_numel_" + id;

    if not(exist(save_dir, 'dir'))
        mkdir(save_dir);
    end

    FONT = 'Times';
    fs = 8;

    % 全体の解のデータ
    data_dir3c = "GS_input/z2033_r602/";
    vars3c = load(data_dir3c + "vars");
    env3c = vars3c.param;
    flux3c = vars3c.psi0;
    psi_v_3c = flux3c';

    % 容器の大きさに拡張したい部分的な解のデータ
    data_dir = "GS_input/" + dirname + "/";
    vars = load(data_dir + "vars");
    env = vars.param;
    flux = vars.psi;
    flux0 = vars.psi0;
    psiorig0 = flux';
    psi0 = flux' - flux0';

    % グラフの領域2020/12/21
    z = linspace(env3c.zmin, env3c.zmax, env3c.Nz);
    r = linspace(env3c.rmin, env3c.rmax, env3c.Nr);

    if env.Nz == 2033
        % 拡張しなくて良い

        psi = flux';
        figure()
        subplot(1, 5, 1)
        v = linspace(-20, 20, 21);
        contour(r, z, flipud(psi' * 1000), v, 'r')
        psi222 = 0;
        jt = 0;
    else
        % 拡張する

        % 片方に数値解をつくった磁場のグラフ2020/12/21
        A = size(psi0);
        psi000 = zeros(env3c.Nr, env3c.Nz - A(1, 2));
        psi00 = [psi0, psi000];
        figure
        ah = subplot(1, 5, 1);
        v = linspace(-20, 20, 21);
        psiorig0 = [psiorig0, psi000];
        contour(r, z, flipud(psiorig0' * 1000), v, 'r')
        ah = subplot(1, 5, 2);
        contour(r, z, flipud(psi00' * 1000), v, 'r')
        % 片方に数値解をつくった磁場のグラフここまで2020/12/21

        psi = psi_v_3c * 0;
        psi(1:env.Nr, 1:env.Nz) = psi0;
        psi(1:env.Nr, end - env.Nz + 1:end) = psi(1:env.Nr, end - env.Nz + 1:end) + psi0(:, end:-1:1);
        % 反対側にコピーした磁場のグラフ2020/12/21
        ah = subplot(1, 5, 3);
        contour(r, z, psi' * 1000, v, 'r')
        % 反対側にコピーした磁場のグラフここまで2020/12/21
        psi = psi + psi_v_3c;

        %% Calc virtual eddy current on center
        flux = vars.mu_jt;
        jt = psi_v_3c * 0;
        jt(1:env.Nr, 1:env.Nz) = jt(1:env.Nr, 1:env.Nz) + flux(1:env.Nz, 1:env.Nr)' / (4 * pi * 1e-7);
        jt(1:env.Nr, end:-1:end - env.Nz + 1) = jt(1:env.Nr, end:-1:end - env.Nz + 1) + ...
            flux(1:env.Nz, 1:env.Nr)' / (4 * pi * 1e-7);
        delz = z(2) - z(1);
        delr = r(2) - r(1);
        jt_center = -squeeze(psi0(2:env.Nr - 1, env.Nz - 1)) / delz ./ r(2:env.Nr - 1)' / 2 / pi / 4 / pi / 1e-7 * delr;
        [rr zz] = meshgrid(r, z(1:ceil(length(z) / 2)));
        psi_virtualj = rr * 0;

        const = 2 * pi * 4 * pi * 1e-7;

        for k = 1:length(jt_center)
            kk1 = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz - z(env.Nz - 1)).^2);
            [K1, E1] = ellipke(kk1);
            kk2 = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz + z(env.Nz - 1)).^2);
            [K2, E2] = ellipke(kk2);
            denominator = ((1 ./ (pi * sqrt(kk1)) .* ((1 - kk1 / 2) .* K1 - E1)) + (1 ./ (pi * sqrt(kk2)) .* ((1 - kk2 / 2) .* K2 - E2)));
            psi_virtualj = psi_virtualj + const * jt_center(k) .* sqrt(rr * r(k + 1)) .* denominator;
        end

        psi222 = psi_v_3c * 0;
        psi222(1:env3c.Nr, 1:(env3c.Nz + 1) / 2) = psi_virtualj(1:(env3c.Nz + 1) / 2, 1:env3c.Nr)';
        psi222(1:env3c.Nr, end:-1:end - (env3c.Nz + 1) / 2 + 1 + 1) = psi_virtualj(1:(env3c.Nz + 1) / 2 - 1, 1:env3c.Nr)';

        % 仮想壁に流れる電流を差し引く
        % 平衡解2020/12/21
        ah = subplot(1, 5, 4);
        contour(r, z, psi' * 1000, v, 'r')
        psi = psi - psi222;
        ah = subplot(1, 5, 5);
        contour(r, z, psi' * 1000, v, 'r')
    end

    %% Calc Bz and Br from psi
    zdiff = z(2:end - 1);
    rdiff = r(2:end - 1);
    Bz = zeros(length(r), length(z));
    Br = zeros(length(r), length(z));

    for j = 1:length(z)
        Bz(2:end - 1, j) = (psi(3:end, j) - psi(1:end - 2, j)) ./ (r(3:end) - r(1:end - 2))' / 2 / pi ./ r(2:end - 1)';
    end

    for j = 1:length(r)
        Br(j, 2:end - 1) = -(psi(j, 3:end) - psi(j, 1:end - 2)) ./ (z(3:end) - z(1:end - 2)) / 2 / pi / r(j);
    end

    %% CCSdata
    datanum = 0;

    % 602*2033
    % 1017 ... Nzの中央
    % 1016 ...  Nzの中央-1
    % 2028 ... Nzから5点内側
    % 6    ... Nzから５点内側

    % 2021/04/28 UTST装置のセンサー配置
    % save("psiGSplotTest", "psi");

    for i = 1:length(z(1017:2028))
        r_CCS(i) = r(6);
        z_CCS(i) = z(1016 + i);
        psi_CCS(i) = psi(6, 1016 + i);
        Bz_CCS(i) = Bz(6, 1016 + i);
        Br_CCS(i) = Br(6, 1016 + i);
    end

    datanum = length(r_CCS);

    for i = 1:length(r(7:499))
        r_CCS(datanum + i) = r(6 + i);
        z_CCS(datanum + i) = z(2028);
        psi_CCS(datanum + i) = psi(6 + i, 2028);
        Bz_CCS(datanum + i) = Bz(6 + i, 2028);
        Br_CCS(datanum + i) = Br(6 + i, 2028);
    end

    datanum = length(r_CCS);

    for i = 1:length(z(1302:2027))
        r_CCS(datanum + i) = r(499);
        z_CCS(datanum + i) = z(2028 - i);
        psi_CCS(datanum + i) = psi(499, 2028 - i);
        Bz_CCS(datanum + i) = Bz(499, 2028 - i);
        Br_CCS(datanum + i) = Br(499, 2028 - i);
    end

    datanum = length(r_CCS);

    for i = 1:length(r(500:597))
        r_CCS(datanum + i) = r(499 + i);
        z_CCS(datanum + i) = z(1302);
        psi_CCS(datanum + i) = psi(499 + i, 1302);
        Bz_CCS(datanum + i) = Bz(499 + i, 1302);
        Br_CCS(datanum + i) = Br(499 + i, 1302);
    end

    datanum = length(r_CCS);

    for i = 1:length(z(1017:1301))
        r_CCS(datanum + i) = r(597);
        z_CCS(datanum + i) = z(1302 - i);
        psi_CCS(datanum + i) = psi(597, 1302 - i);
        Bz_CCS(datanum + i) = Bz(597, 1302 - i);
        Br_CCS(datanum + i) = Br(597, 1302 - i);
    end

    if saveflag == "origin"
        % テキストファイルで保存するならこっち
        % matファイルで保存するなら最後

        if realFlag
            makeRealSenposCCSdata(psi, Bz, Br, save_dir);
            % error('test');
            % ここまで
        else
            datanum = length(r_CCS);

            fp = fopen(save_dir + '/Sensor_.txt', 'w');
            fprintf(fp, 'r[m]\tz[m]\tpsi[Wb]\tBz[T]\tBr[T]\n');

            for i = 1:30:datanum
                fprintf(fp, '%f\t%f\t%f\t%f\t%f\n', r_CCS(i), z_CCS(i), psi_CCS(i), Bz_CCS(i), Br_CCS(i));
            end

            fclose(fp);
        end

        fp = fopen(save_dir + '/CCS.txt', 'w');
        fprintf(fp, 'r[m]\tz[m]\tpsi[Wb]\tBz[T]\tBr[T]\n');

        for i = 1:datanum
            fprintf(fp, '%f\t%f\t%f\t%f\t%f\n', r_CCS(i), z_CCS(i), psi_CCS(i), Bz_CCS(i), Br_CCS(i));
        end

        fclose(fp);

        fp = fopen(save_dir + '/FluxProfile_2D.txt', 'w');

        for i = 1:length(z)

            for j = 1:length(r) - 1
                fprintf(fp, '%f\t', psi(j, i));
            end

            fprintf(fp, '%f\n', psi(length(r), i));
        end

        fclose(fp);
    end

    if dispFigure

        fh = figure;
        set(fh, 'position', [50 50 600 600], 'Name', ['psi contour'], 'NumberTitle', 'off');

        ah = subplot(3, 1, 1);
        set(ah, 'box', 'on', 'FontSize', fs, 'FontName', FONT);

        hold on
        pcolor(z, r, sqrt(Bz.^2 + Br.^2) * 1e3)
        axis equal
        shading interp

        ca = [-100 100];
        caxis(ca);

        v = linspace(-100, 100, 31);
        contour(z, r, sqrt(Bz.^2 + Br.^2) * 10000, v, 'k')
        v = linspace(0, 5, 11);
        contour(z, r, sqrt(Bz.^2 + Br.^2) * 10000, v, 'k')
        plot([env3c.zmin env3c.zmin env3c.z1 env3c.z2 env3c.z3 env3c.z4 env3c.zmax env3c.zmax], [env3c.rmin env3c.rmirror env3c.rmirror env3c.rmax env3c.rmax env3c.rmirror env3c.rmirror env3c.rmin], 'r', 'LineWidth', 2);
        title('|Bp| [mT]');
        set(ah, 'Box', 'on', 'FontSize', fs, 'FontName', FONT, 'FontWeight', 'bold'); %,'xlim', [-1.2 1.2], 'ylim', [0 0.8]);

        ah = subplot(3, 1, 2);
        set(ah, 'box', 'on', 'FontSize', fs, 'FontName', FONT);

        hold on
        pcolor(z, r, psi * -1000)
        axis equal
        shading interp

        ca = [-5 5];
        caxis(ca);

        v = linspace(-15, 15, 11);
        contour(z, r, psi * -1000, v, 'k')
        v = linspace(-5, 5, 11);
        contour(z, r, psi * -1000, v, 'w')
        %contour(z,r,psi*-1000,-2:0.1:0,'r')
        plot([env3c.zmin env3c.zmin env3c.z1 env3c.z2 env3c.z3 env3c.z4 env3c.zmax env3c.zmax], [env3c.rmin env3c.rmirror env3c.rmirror env3c.rmax env3c.rmax env3c.rmirror env3c.rmirror env3c.rmin], 'r', 'LineWidth', 2);

        % plot(z_PS,r_PS,'r','LineWidth',2);
        % plot(z_LCFS,r_LCFS,'g','LineWidth',2);
        title('poloidal flux [mWb]');
        set(ah, 'Box', 'on', 'FontSize', fs, 'FontName', FONT, 'FontWeight', 'bold'); %,'xlim', [-1.2 1.2], 'ylim', [0 0.8]);

        ah = subplot(3, 1, 3);
        set(ah, 'box', 'on', 'FontSize', fs, 'FontName', FONT);
    end

    if dispFigure
        hold on
        pcolor(z, r, -jt / 1e3)
        axis equal
        shading interp

        ca = [-500 500];
        caxis(ca);

        v = linspace(-500, 500, 11);
        contour(z, r, -jt / 1e3, v, 'k')
        v = linspace(-5, 5, 11);
        contour(z, r, -jt / 1e3, v, 'w')
        contour(z, r, -jt / 1e3, 10:1:-10, 'r')
        plot([env3c.zmin env3c.zmin env3c.z1 env3c.z2 env3c.z3 env3c.z4 env3c.zmax env3c.zmax], [env3c.rmin env3c.rmirror env3c.rmirror env3c.rmax env3c.rmax env3c.rmirror env3c.rmirror env3c.rmin], 'r', 'LineWidth', 2);
        title('current density [kA/m2]');
        set(ah, 'Box', 'on', 'FontSize', fs, 'FontName', FONT, 'FontWeight', 'bold'); %,'xlim', [-1.2 1.2], 'ylim', [0 0.8]);
    end

    %% Calc eddy current
    psi_eddy = psi - psi_v_3c + psi222;

    delz = z(2) - z(1);
    delr = r(2) - r(1);

    jt1 = (-psi_eddy(601, 1017:1307) / delr / r(602) ...
        - psi_eddy(601, 1017:1307) / 2 / r(602)^2) ...
        / 2 / pi / 4 / pi / 1e-7/1e6;

    jt2 = -squeeze(psi_eddy(600:-1:502, 1306)) / delz ./ r(600:-1:502)' ...
        / 2 / pi / 4 / pi / 1e-7/1e6;

    jt3 = (-psi_eddy(502, 1306:2032) / delr / r(503) ...
        - psi_eddy(502, 1306:2032) / 2 / r(503)^2) ...
        / 2 / pi / 4 / pi / 1e-7/1e6;

    jt4 = -squeeze(psi_eddy(501:-1:2, 2032)) / delz ./ r(501:-1:2)' ...
        / 2 / pi / 4 / pi / 1e-7/1e6;

    jt5 = (-psi_eddy(2, 2031:-1:1017) / delr / r(1) ... % d^2 psi/dr^2
        + psi_eddy(2, 2031:-1:1017) / 2 / r(1)^2) ... % d psi/dr
        / 2 / pi / 4 / pi / 1e-7/1e6;

    L1 = delz * (1:length(jt1));
    L2 = delz * length(jt1) + delr * (1:length(jt2));
    L3 = delz * length(jt1) + delr * length(jt2) + delz * (1:length(jt3));
    L4 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * (1:length(jt4));
    L5 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * (1:length(jt5));
    L6 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
        delz * (1:length(jt5));
    L7 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
        delz * length(jt5) + delr * (1:length(jt4));
    L8 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
        delz * length(jt5) + delr * length(jt4) + delz * (1:length(jt3));
    L9 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
        delz * length(jt5) + delr * length(jt4) + delz * length(jt3) + delr * (1:length(jt2));
    L10 = delz * length(jt1) + delr * length(jt2) + delz * length(jt3) + delr * length(jt4) + delz * length(jt5) + ...
        delz * length(jt5) + delr * length(jt4) + delz * length(jt3) + delr * length(jt2) + delz * (1:length(jt1));

    Length = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10];
    jeddy = [jt1 jt2' jt3 jt4' jt5 jt5(end:-1:1) jt4(end:-1:1)' jt3(end:-1:1) jt2(end:-1:1)' jt1(end:-1:1)];

    if dispFigure
        fh = figure;
        hold on
        plot(L1, jt1, 'r');
        plot(L2, jt2, 'r');
        plot(L3, jt3, 'r');
        plot(L4, jt4, 'r');
        plot(L5, jt5, 'r');
        plot(L6, jt5(end:-1:1), 'r');
        plot(L7, jt4(end:-1:1), 'r');
        plot(L8, jt3(end:-1:1), 'r');
        plot(L9, jt2(end:-1:1), 'r');
        plot(L10, jt1(end:-1:1), 'r');
        plot(Length, cumtrapz(Length, jeddy))
        xlabel('Distance (m)');
        ylabel('Eddy current density (MA/m)');
    end

    if saveflag == "origin"
        fp = fopen(save_dir + '/jeddy.txt', 'w');

        for j = 1:length(Length)
            fprintf(fp, '%f\t%f\n', Length(j), jeddy(j));
        end

        fclose(fp);
    elseif saveflag == "modan"
        r_CCS = [r_CCS fliplr(r_CCS)];
        z_CCS = [z_CCS -fliplr(z_CCS)];
        Br_CCS = [Br_CCS -fliplr(Br_CCS)];
        Bz_CCS = [Bz_CCS fliplr(Bz_CCS)];
        psi_CCS = [psi_CCS fliplr(psi_CCS)];
        save([save_dir + '/merged.mat'], "r_CCS", "z_CCS", "psi_CCS", "Bz_CCS", "Br_CCS", "jeddy", "psi");
    end

end
