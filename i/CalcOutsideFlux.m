% suffix = "test"
% save_dir = "GSPLOT_OUTPUT";

% fp = fopen(save_dir+'/jeddy_'+suffix+'.txt', 'w')
% fprintf(fp, 'a');

% SEN0 = dlmread('i/sensorCoordinate0.txt');
% r = SEN0(:, 1);
% z = SEN0(:, 2);
% figure()
% plot(r, z, 'o')

% 容器の外側の磁場を計算
datafile_dir = "INPUT/z2033_r602/";

param.file0 = datafile_dir + "def_mesh.dat";
param.file1 = datafile_dir + "def_slow.dat";

% ファイルから値の読み込み
A = load(param.file0);
param.Nz = A(1, 1);
param.Nr = A(2, 1);
param.zmin = A(3, 1);
param.zmax = A(4, 1);
param.rmin = A(5, 1);
param.rmax = A(6, 1);
param.z1 = A(7, 1);
param.z2 = A(8, 1);
param.z3 = A(9, 1);
param.z4 = A(10, 1);
param.Nr = 800;
param.rmax = 0.8880;

% コイルに関するデータの読み込み
f = load(param.file1);
param.coil_turn = f(1, :);
param.coil_Ic = f(2, :) .* f(1, :);
param.coil_z = f(3, :);
param.coil_r = f(4, :);
param.ncoil = size(param.coil_turn, 2);

psi = zeros(param.Nz, param.Nr);
[param, psi] = cal_psi(param, psi);

z = linspace(param.zmin, param.zmax, param.Nz);
r = linspace(param.rmin, param.rmax, param.Nr);
v = linspace(-50, 50, 21);

figure()

contour(r, z, psi * 1000, v)
% save("psi800", "psi")

function [param, psi] = cal_psi(param, psi)

    Mu = 4 * pi * 1.0e-7; % 真空の透磁率
    mtrxz = (0:1:param.Nz - 1)' .* param.delz + param.zmin;
    mtrxr = (0:1:param.Nr - 1) .* param.delr + param.rmin;

    for k = 1:param.ncoil
        % (z*rの大きさの配列を作成)
        tmp1 = (param.coil_z(1, k) - mtrxz).^2 * ones(1, param.Nr) + ones(param.Nz, 1) * (param.coil_r(1, k) + mtrxr).^2;
        kk = (4 .* ones(param.Nz, 1) * mtrxr .* param.coil_r(1, k)) ./ tmp1;
        [K, E] = ellipke(kk);
        psi = psi + 2 * pi * Mu * param.coil_Ic(1, k) ./ (pi * (kk).^0.5) .* sqrt(param.coil_r(1, k) .* mtrxr) .* ...
            ((1 - kk / 2) .* K - E);
    end

end
