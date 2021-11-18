dirname = 'z1000_r602';
dirname = 'z2033_r602';

FONT = 'Times';
fs = 8;

% 全体の解のデータ
data_dir3c = "GS_input/z2033_r602/";
% data_dir3c = "GS_input/z201_r101/";
vars3c = load(data_dir3c + "vars");
env3c = vars3c.param;
flux3c = vars3c.psi0;
psi_v_3c = flux3c';

% 容器の大きさに拡張したい部分的な解のデータ
data_dir = "GS_input/" + dirname + "/";
vars = load(data_dir + "vars");
env = vars.param;
flux = vars.psi; % GSコードで求めた磁束
flux0 = vars.psi0; % EFコイルが作る磁束
psi0 = flux' - flux0'; % コイルが作る磁束を差し引いた磁束
psiorig0 = flux';

% グラフの領域2020/12/21
z = linspace(env3c.zmin, env3c.zmax, env3c.Nz);
r = linspace(env3c.rmin, env3c.rmax, env3c.Nr);

delz = z(2) - z(1);
delr = r(2) - r(1);

if env.Nz == env3c.Nz
    % 拡張しなくて良い

    psi = flux';
    figure()
    subplot(1, 6, 1)
    v = linspace(-20, 20, 21);
    contour(r, z, flipud(psi' * 1000), v, 'r')
    psi222 = 0;
    jt = 0;
    psi_nocoil = psi0;
else
    % 拡張する

    % 片方に数値解をつくった磁場のグラフ2020/12/21
    A = size(psi0);
    psi000 = zeros(env3c.Nr, env3c.Nz - A(1, 2));
    psi00 = [psi0, psi000];
    figure
    ah = subplot(1, 6, 1);
    v = linspace(-20, 20, 21);
    psiorig0 = [psiorig0, psi000];
    contour(r, z, flipud(psiorig0' * 1000), v, 'r')
    ah = subplot(1, 6, 2);
    contour(r, z, flipud(psi00' * 1000), v, 'r')
    % 片方に数値解をつくった磁場のグラフここまで2020/12/21

    psi = psi_v_3c * 0;
    psi(1:env.Nr, 1:env.Nz) = psi0;
    psi(1:env.Nr, end - env.Nz + 1:end) = psi(1:env.Nr, end - env.Nz + 1:end) + psi0(:, end:-1:1);
    % 反対側にコピーした磁場のグラフ2020/12/21
    ah = subplot(1, 6, 3);
    contour(r, z, psi' * 1000, v, 'r')
    % 外部コイルが作る磁束を足した磁束
    psi_nocoil = psi;
    psi = psi + psi_v_3c;
    ah = subplot(1, 6, 4);
    contour(r, z, psi' * 1000, v, 'r')
    %

    % error('error description', A1)
    %% Calc virtual eddy current on center
    flux = vars.mu_jt;
    jt = psi_v_3c * 0;
    jt(1:env.Nr, 1:env.Nz) = jt(1:env.Nr, 1:env.Nz) + flux(1:env.Nz, 1:env.Nr)' / (4 * pi * 1e-7);
    jt(1:env.Nr, end:-1:end - env.Nz + 1) = jt(1:env.Nr, end:-1:end - env.Nz + 1) + ...
        flux(1:env.Nz, 1:env.Nr)' / (4 * pi * 1e-7);
    jt_center = -squeeze(psi0(2:env.Nr - 1, env.Nz - 1)) / delz ./ r(2:env.Nr - 1)' / 2 / pi / 4 / pi / 1e-7 * delr;
    size(jt_center)
    [rr zz] = meshgrid(r, z(1:ceil(length(z) / 2)));
    psi_virtualj = rr * 0;

    const = 2 * pi * 4 * pi * 1e-7;

    % error('error description', A1)

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
    psi = psi - psi222;
    ah = subplot(1, 6, 5);
    contour(r, z, psi' * 1000, v, 'r')
end

%% UTST凸部の渦電流を差し引く
% UTST凸部のZ座標（zUが上zLが下）
[M zU_mesh] = min(abs(z - env3c.z3));
[M zL_mesh] = min(abs(z - env3c.z2));
% 仮想壁部の１つ内側の電流密度
jt_wall = -squeeze(psi_nocoil(env.Nr - 1, zL_mesh + 1:zU_mesh - 1)) / delz / r(env.Nr - 1) / 2 / pi / 4 / pi / 1e-7 * delr;
[rr zz] = meshgrid(r, z(1:ceil(length(z) / 2)));
psi_virtualj = rr * 0;
const = 2 * pi * 4 * pi * 1e-7;
i = 1;

for k = zL_mesh + 1:zU_mesh - 1
    kk1 = 4 * rr * r(env.Nr - 1) ./ ((rr + r(env.Nr - 1) + 1e-6).^2 + (zz - z(k)).^2);
    [K1, E1] = ellipke(kk1);
    kk2 = 4 * rr * r(env.Nr - 1) ./ ((rr + r(env.Nr - 1) + 1e-6).^2 + (zz + z(k)).^2);
    [K2, E2] = ellipke(kk2);
    denominator = ((1 ./ (pi * sqrt(kk1)) .* ((1 - kk1 / 2) .* K1 - E1)) + (1 ./ (pi * sqrt(kk2)) .* ((1 - kk2 / 2) .* K2 - E2)));
    psi_virtualj = psi_virtualj + const * jt_wall(i) .* sqrt(rr * r(env.Nr - 1)) .* denominator;
    i = i + 1;
end

psi2222 = psi_v_3c * 0;
psi2222(1:env3c.Nr, 1:(env3c.Nz + 1) / 2) = psi_virtualj(1:(env3c.Nz + 1) / 2, 1:env3c.Nr)';
psi2222(1:env3c.Nr, end:-1:end - (env3c.Nz + 1) / 2 + 1 + 1) = psi_virtualj(1:(env3c.Nz + 1) / 2 - 1, 1:env3c.Nr)';
psi = psi - psi2222;
ah = subplot(1, 6, 6);
contour(r, z, psi' * 1000, v, 'r')

figure()
subplot(1, 3, 1);
contour(r, z, flux * 1000, v)
subplot(1, 3, 2);
contour(r, z, psi' * 1000, v)
subplot(1, 3, 3);
contour(r, z, psi2222' * 1000, v)
%% UTST凸部の渦電流を差し引く計算ここまで
