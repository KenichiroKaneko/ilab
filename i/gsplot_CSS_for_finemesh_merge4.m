
% CCS03cはUTST装置全体の解
% CCS05は（途中に仮想壁を置いた）部分的な解
% 部分的な解を全体に拡張する際に、CCS03cのメッシュデータが必要
% fluxは平衡解、flux0はpsiの初期状態

% data_dirs = ["z2033_r602" ...
% data_dirs = ["z1000_r602" "z0960_r602" "z0920_r602" "z0880_r602" ...
%             "z0840_r602" "z0800_r602" "z0760_r602" "z0720_r602" ...
%             "z0680_r602" "z0640_r602" "z0600_r602" "z0560_r602" ...
%             "z0520_r602" "z0480_r602" "z0440_r602" "z0400_r602" ...
%             "z0360_r602" "z0320_r602" "z0280_r602" "z0240_r602" ...
%             "z0200_r602" "z0160_r602" "z0120_r602" "z0080_r602"];
% data_dirs = ["z0200_r602" "z0160_r602" "z0120_r602" "z0080_r602"];
data_dirs = ["z1000_r602"];
saveflag = 0;
FONT = 'Times';
fs = 8;


figure

for i = 1:length(data_dirs)
    gsplot_CCS(i, data_dirs(i), length(data_dirs));
end 
    
function gsplot_CCS(i, dir_name, data_num)
    % 全体の解のデータ
    data_dir3c = "outputfiles_03c/";
    vars3c = load(data_dir3c + "vars");
    env3c = vars3c.param;
    flux3c = vars3c.psi0;
    psi_v_3c = flux3c';
    % グラフの領域2020/12/21
    z = linspace(env3c.zmin, env3c.zmax, env3c.Nz);
    r = linspace(env3c.rmin, env3c.rmax, env3c.Nr);
    % 部分的な解のデータ
    data_dir = "OUTPUT/" + dir_name;
    vars = load(data_dir + "/vars");
    env = vars.param;
    flux = vars.psi;
    flux0 = vars.psi0;
    psiorig0 = flux';
    psi0 = flux' - flux0';



    % 片方に数値解をつくった磁場のグラフ2020/12/21
    A = size(psi0);
    psi000 = zeros(env3c.Nr, env3c.Nz - A(1, 2));
    psi00 = [psi0, psi000];
    % ah = subplot(1, 5, 1);
    v = linspace(-20, 20, 21);
    psiorig0 = [psiorig0, psi000];
    % contour(r, z, flipud(psiorig0' * 1000), v, 'r')
    % ah = subplot(1, 5, 2);
    % contour(r, z, flipud(psi00' * 1000), v, 'r')
    % 片方に数値解をつくった磁場のグラフここまで2020/12/21

    psi = psi_v_3c * 0;
    psi(1:env.Nr, 1:env.Nz) = psi0;
    psi(1:env.Nr, end - env.Nz + 1:end) = psi(1:env.Nr, end - env.Nz + 1:end) + psi0(:, end:-1:1);
    % 反対側にコピーした磁場のグラフ2020/12/21
    % ah = subplot(1, 5, 3);
    % contour(r, z, psi' * 1000, v, 'r')
    % 反対側にコピーした磁場のグラフここまで2020/12/21
    psi = psi + psi_v_3c;

    flux = vars.mu_jt;
    jt = psi_v_3c * 0;
    jt(1:env.Nr, 1:env.Nz) = jt(1:env.Nr, 1:env.Nz) + flux(1:env.Nz, 1:env.Nr)' / (4 * pi * 1e-7);
    jt(1:env.Nr, end:-1:end - env.Nz + 1) = jt(1:env.Nr, end:-1:end - env.Nz + 1) + ...
        flux(1:env.Nz, 1:env.Nr)' / (4 * pi * 1e-7);

    %% Calc virtual eddy current on center
    delz = z(2) - z(1);
    delr = r(2) - r(1);

    jt_center = -squeeze(psi0(2:env.Nr - 1, env.Nz - 1)) / delz ./ r(2:env.Nr - 1)' / 2 / pi / 4 / pi / 1e-7 * delr;

    [rr zz] = meshgrid(r, z(1:ceil(length(z) / 2)));

    psi_virtualj = rr * 0;

    for k = 1:length(jt_center)
        kk = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz - z(env.Nz - 1)).^2);
        [K, E] = ellipke(kk);
        %    K(~isfinite(K))=0;
        %    E(~isfinite(E))=0;
        psi_virtualj = psi_virtualj + 2 * pi * 4 * pi * 1e-7 * jt_center(k) ./ (pi * sqrt(kk)) .* sqrt(rr * r(k + 1)) .* ((1 - kk / 2) .* K - E);
        kk = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz + z(env.Nz - 1)).^2);
        [K, E] = ellipke(kk);
        %    K(~isfinite(K))=0;
        %    E(~isfinite(E))=0;
        psi_virtualj = psi_virtualj + 2 * pi * 4 * pi * 1e-7 * jt_center(k) ./ (pi * sqrt(kk)) .* sqrt(rr * r(k + 1)) .* ((1 - kk / 2) .* K - E);
    end

    psi222 = psi_v_3c * 0;
    psi222(1:env3c.Nr, 1:(env3c.Nz + 1) / 2) = psi_virtualj(1:(env3c.Nz + 1) / 2, 1:env3c.Nr)';
    psi222(1:env3c.Nr, end:-1:end - (env3c.Nz + 1) / 2 + 1 + 1) = psi_virtualj(1:(env3c.Nz + 1) / 2 - 1, 1:env3c.Nr)';
    % 平衡解2020/12/21
    % contour(r, z, psi' * 1000, v, 'r')
    psi = psi - psi222;

    ah = subplot(3, 8, i);
    contour(r, z, psi' * 1000, v, 'r')
end