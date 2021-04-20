function gsplot_CCS_for_finemesh_merge(suffix)
    %clear  

    saveflag = 1;

    FONT = 'Times';
    fs = 8;
    dir = (['CCS03c_inputs/']);
    % dir = (['outputfiles/']);

    env = load([dir 'equiv.txt']);
    flux = load([dir 'flux0.txt']);
    psi_v_3c = flux';   

    %suffix='9cEF150';
    dir = (['CCS05c_inputs/']);
    %dir=(['c:/cygwin64/home/inomoto/prog/GS_for_CCS/UTST_CCS04c_double_finemesh/outputfiles/']);

    env_tmp = load([dir 'equiv.txt']);

    flux = load([dir 'flux.txt']);
    psi0 = flux';
    psiorig0 = psi0;
    flux = load([dir 'flux0.txt']);
    psi0 = psi0 - flux';


    % グラフの領域2020/12/21
    z = linspace(env(3), env(4), env(1));
    r = linspace(env(5), env(6), env(2));

    % 片方で数値解をつくるやつ2020/12/21
    A = size(psi0);
    psi000 = zeros(602, 2033 - A(1, 2));
    psi00 = [psi0, psi000];
    figure
    ah = subplot(1,5,1);
    v = linspace(-20, 20, 21);
    psiorig0 = [psiorig0, psi000];
    contour(r, z, flipud(psiorig0' * 1000), v, 'r')
    ah = subplot(1,5,2);
    contour(r, z, flipud(psi00' * 1000), v, 'r')

    %psi1=psi0(:,end-1:-1:1);
    %psi=[psi0 psi1]+psi_v_3c;
    psi = psi_v_3c * 0;
    psi(1:env_tmp(2), 1:env_tmp(1)) = psi0;
    psi(1:env_tmp(2), end - env_tmp(1) + 1:end) = psi(1:env_tmp(2), end - env_tmp(1) + 1:end) + psi0(:, end:-1:1);
    % 反対側にコピー2020/12/21
    ah = subplot(1,5,3);
    contour(r, z, psi' * 1000, v, 'r')
    %psi(1:602,1521+1:end-1521)=psi_v_3c(1:602,1521+1:end-1521);
    psi = psi + psi_v_3c;


    flux = load([dir 'mujt.txt']);
    jt = psi_v_3c * 0;
    jt(1:env_tmp(2), 1:env_tmp(1)) = jt(1:env_tmp(2), 1:env_tmp(1)) + flux(1:env_tmp(1), 1:env_tmp(2))' / (4 * pi * 1e-7);
    jt(1:env_tmp(2), end:-1:end - env_tmp(1) + 1) = jt(1:env_tmp(2), end:-1:end - env_tmp(1) + 1) + flux(1:env_tmp(1), 1:env_tmp(2))' / (4 * pi * 1e-7);


    %% Calc virtual eddy current on center
    delz = z(2) - z(1);
    delr = r(2) - r(1);

    jt_center = -squeeze(psi0(2:env_tmp(2) - 1, env_tmp(1) - 1)) / delz ./ r(2:env_tmp(2) - 1)' / 2 / pi / 4 / pi / 1e-7 * delr;

    [rr zz] = meshgrid(r, z(1:ceil(length(z) / 2)));

    psi_virtualj = rr * 0;

    for k = 1:length(jt_center)
        kk = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz - z(env_tmp(1) - 1)).^2);
        [K, E] = ellipke(kk);
        %    K(~isfinite(K))=0;
        %    E(~isfinite(E))=0;
        psi_virtualj = psi_virtualj + 2 * pi * 4 * pi * 1e-7 * jt_center(k) ./ (pi * sqrt(kk)) .* sqrt(rr * r(k + 1)) .* ((1 - kk / 2) .* K - E);
        kk = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz + z(env_tmp(1) - 1)).^2);
        [K, E] = ellipke(kk);
        %    K(~isfinite(K))=0;
        %    E(~isfinite(E))=0;
        psi_virtualj = psi_virtualj + 2 * pi * 4 * pi * 1e-7 * jt_center(k) ./ (pi * sqrt(kk)) .* sqrt(rr * r(k + 1)) .* ((1 - kk / 2) .* K - E);
    end

    psi222 = psi_v_3c * 0;
    psi222(1:env(2), 1:(env(1) + 1) / 2) = psi_virtualj(1:(env(1) + 1) / 2, 1:env(2))';
    psi222(1:env(2), end:-1:end - (env(1) + 1) / 2 + 1 + 1) = psi_virtualj(1:(env(1) + 1) / 2 - 1, 1:env(2))';
    % 平衡解2020/12/21
    ah = subplot(1,5,4);
    contour(r, z, psi' * 1000, v, 'r')
    psi = psi - psi222;
    ah = subplot(1, 5, 5);
    contour(r, z, psi' * 1000, v, 'r')
    % error('fin')

    % jt_center=jt_center*0-15;
    % psi_virtualj=rr*0;
    % for k=1:length(jt_center)
    %     kk = 4 * rr*r(k+1)./((rr+r(k+1)+1e-6).^2+(zz-0).^2);
    %     [K,E]=ellipke(kk);
    % %    K(~isfinite(K))=0;
    % %    E(~isfinite(E))=0;
    %     psi_virtualj= psi_virtualj+2*pi*4*pi*1e-7*jt_center(k)./(pi*sqrt(kk)).*sqrt(rr*r(k+1)).*((1-kk/2).*K-E);
    % end
    % psi222=psi_v_3c*0;
    % psi222(1:602,1:1017)=psi_virtualj(1:1017,1:602)';
    % psi222(1:602,end:-1:end-1017+1+1)=psi_virtualj(1:1017-1,1:602)';
    % psi=psi-psi222;

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

    datanum = length(r_CCS);
    % figure()
    % plot( r_CCS, z_CCS, 'o');
    % error('fin')

    if (saveflag)
        fp = fopen(['CCS_' suffix '.txt'], 'w');
        fprintf(fp, 'r[m]\tz[m]\tpsi[Wb]\tBz[T]\tBr[T]\n');

        for i = 1:datanum
            fprintf(fp, '%f\t%f\t%f\t%f\t%f\n', r_CCS(i), z_CCS(i), psi_CCS(i), Bz_CCS(i), Br_CCS(i));
        end

        fclose(fp);

        fp = fopen(['CCS_thin_' suffix '.txt'], 'w');
        fprintf(fp, 'r[m]\tz[m]\tpsi[Wb]\tBz[T]\tBr[T]\n');

        for i = 1:5:datanum
            fprintf(fp, '%f\t%f\t%f\t%f\t%f\n', r_CCS(i), z_CCS(i), psi_CCS(i), Bz_CCS(i), Br_CCS(i));
        end

        fclose(fp);

        fp = fopen(['flux_' suffix '.txt'], 'w');
        %fprintf(fp,'r[m]\tz[m]\tpsi[Wb]\tBz[T]\tBr[T]\n');
        for i = 1:length(z)

            for j = 1:length(r) - 1
                fprintf(fp, '%f\t', psi(j, i));
            end

            fprintf(fp, '%f\n', psi(length(r), i));
        end

        fclose(fp);
    end

    % %% Plasma surface & LCFS
    %
    % num_PS1=1;
    % num_LCFS1=1;
    % for j=1:length(z)
    %     if(min(psi(:,j))<0)
    %         i_r_psimin=find(psi(:,j)==min(psi(:,j)));
    %         z_PS(num_PS1)=z(j);
    %         r_PS(num_PS1)=interp1(psi(1:i_r_psimin(1),j),r(1:i_r_psimin(1)),0);
    %         num_PS1=num_PS1+1;
    %     end
    %     if(min(psi(:,j))<5e-3)
    %         i_r_psimin=find(psi(:,j)==min(psi(:,j)));
    %         z_LCFS(num_LCFS1)=z(j);
    %         r_LCFS(num_LCFS1)=interp1(psi(1:i_r_psimin(1),j),r(1:i_r_psimin(1)),5e-3);
    %         num_LCFS1=num_LCFS1+1;
    %     end
    % end
    %
    % num_PS2=num_PS1;
    % num_LCFS2=num_LCFS1;
    % for j=length(z):-1:1
    %     if(min(psi(:,j))<0)
    %         i_r_psimin=find(psi(:,j)==min(psi(:,j)));
    %         z_PS(num_PS2)=z(j);
    %         r_PS(num_PS2)=interp1(psi(i_r_psimin(end):length(r),j),r(i_r_psimin(end):length(r)),0);
    %         num_PS2=num_PS2+1;
    %     end
    %     if(min(psi(:,j))<5e-3)
    %         i_r_psimin=find(psi(:,j)==min(psi(:,j)));
    %         z_LCFS(num_LCFS2)=z(j);
    %         r_LCFS(num_LCFS2)=interp1(psi(i_r_psimin(end):length(r),j),r(i_r_psimin(end):length(r)),5e-3);
    %         num_LCFS2=num_LCFS2+1;
    %     end
    % end
    %
    % num_PS3=num_PS2;
    % num_LCFS3=num_LCFS2;
    % for i=1:length(r)
    %     if((r(i)>r_PS(num_PS1-2))&&(r(i)<r_PS(num_PS1+1)))
    %         r_PS(num_PS3)=r(i);
    %         z_PS(num_PS3)=interp1(psi(i,floor(length(z)/2)+1:length(z)),z(floor(length(z)/2)+1:length(z)),0);
    %         num_PS3=num_PS3+1;
    %     end
    %     if((r(i)>r_LCFS(num_LCFS1-2))&&(r(i)<r_LCFS(num_LCFS1+1)))
    %         r_LCFS(num_LCFS3)=r(i);
    %         z_LCFS(num_LCFS3)=interp1(psi(i,floor(length(z)/2)+1:length(z)),z(floor(length(z)/2)+1:length(z)),5e-3);
    %         num_LCFS3=num_LCFS3+1;
    %     end
    % end
    % num_PS4=num_PS3;
    % num_LCFS4=num_LCFS3;
    % for i=length(r):-1:1
    %     if((r(i)>r_PS(num_PS1-2))&&(r(i)<r_PS(num_PS1+1)))
    %         r_PS(num_PS4)=r(i);
    %         z_PS(num_PS4)=interp1(psi(i,1:floor(length(z)/2)-1),z(1:floor(length(z)/2)-1),0);
    %         num_PS4=num_PS4+1;
    %     end
    %     if((r(i)>r_LCFS(num_LCFS1-2))&&(r(i)<r_LCFS(num_LCFS1+1)))
    %         r_LCFS(num_LCFS4)=r(i);
    %         z_LCFS(num_LCFS4)=interp1(psi(i,1:floor(length(z)/2)-1),z(1:floor(length(z)/2)-1),5e-3);
    %         num_LCFS4=num_LCFS4+1;
    %     end
    % end
    % r_PS=[r_PS(2:num_PS1-2) r_PS(num_PS2:num_PS3-1) r_PS(num_PS1+1:num_PS2-2) r_PS(num_PS3:num_PS4-1)];
    % z_PS=[z_PS(2:num_PS1-2) z_PS(num_PS2:num_PS3-1) z_PS(num_PS1+1:num_PS2-2) z_PS(num_PS3:num_PS4-1)];
    % r_LCFS=[r_LCFS(1:num_LCFS1-1) r_LCFS(num_LCFS2:num_LCFS3-1) r_LCFS(num_LCFS1:num_LCFS2-1) r_LCFS(num_LCFS3:num_LCFS4-1)];
    % z_LCFS=[z_LCFS(1:num_LCFS1-1) z_LCFS(num_LCFS2:num_LCFS3-1) z_LCFS(num_LCFS1:num_LCFS2-1) z_LCFS(num_LCFS3:num_LCFS4-1)];
    %
    % %save('PS','r_PS','z_PS','r_LCFS','z_LCFS');
    % fp=fopen('PS_5c.txt','w');
    % for j=1:length(z_PS)
    %     fprintf(fp,'%f\t%f\n',r_PS(j),z_PS(j));
    % end
    % fclose(fp);
    %
    % fp=fopen('LCFS_5c.txt','w');
    % for j=1:length(z_LCFS)
    %     fprintf(fp,'%f\t%f\n',r_LCFS(j),z_LCFS(j));
    % end
    % fclose(fp);
    %

    %%
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
    plot([env(3) env(3) env(7) env(8) env(9) env(10) env(4) env(4)], [env(5) env(11) env(11) env(6) env(6) env(11) env(11) env(5)], 'r', 'LineWidth', 2);
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
    plot([env(3) env(3) env(7) env(8) env(9) env(10) env(4) env(4)], [env(5) env(11) env(11) env(6) env(6) env(11) env(11) env(5)], 'r', 'LineWidth', 2);

    % plot(z_PS,r_PS,'r','LineWidth',2);
    % plot(z_LCFS,r_LCFS,'g','LineWidth',2);
    title('poloidal flux [mWb]');
    set(ah, 'Box', 'on', 'FontSize', fs, 'FontName', FONT, 'FontWeight', 'bold'); %,'xlim', [-1.2 1.2], 'ylim', [0 0.8]);

    ah = subplot(3, 1, 3);
    set(ah, 'box', 'on', 'FontSize', fs, 'FontName', FONT);

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
    plot([env(3) env(3) env(7) env(8) env(9) env(10) env(4) env(4)], [env(5) env(11) env(11) env(6) env(6) env(11) env(11) env(5)], 'r', 'LineWidth', 2);
    title('current density [kA/m2]');
    set(ah, 'Box', 'on', 'FontSize', fs, 'FontName', FONT, 'FontWeight', 'bold'); %,'xlim', [-1.2 1.2], 'ylim', [0 0.8]);

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

    jt5 = (-psi_eddy(2, 2031:-1:1017) / delr / r(1) ...% d^2 psi/dr^2
    + psi_eddy(2, 2031:-1:1017) / 2 / r(1)^2) ...% d psi/dr
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

    xlabel('Distance (m)');
    ylabel('Eddy current density (MA/m)');

    Length = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10];
    jeddy = [jt1 jt2' jt3 jt4' jt5 jt5(end:-1:1) jt4(end:-1:1)' jt3(end:-1:1) jt2(end:-1:1)' jt1(end:-1:1)];

    if (saveflag)
        fp = fopen(['jeddy_' suffix '.txt'], 'w');

        for j = 1:length(Length)
            fprintf(fp, '%f\t%f\n', Length(j), jeddy(j));
        end

        fclose(fp);
    end

    fh = figure;
    hold on
    plot(Length, cumtrapz(Length, jeddy))

    if (saveflag)
        save(['merged_' suffix '.mat'], 'env', 'delr', 'delz', 'psi', 'psi_v_3c', 'jt', 'r', 'rr', 'z', 'zz');
    end
