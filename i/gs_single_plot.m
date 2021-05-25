% 数値解のみをプロットする

% data_dirs = ["outputfiles_03c/" "outputfiles_04c/" "outputfiles_05c/" "outputfiles_test/"];
% data_dirs = ["outputfiles_test/"];
% data_dirs = ["z2033_r602" ...
    %             "z1000_r602" "z0960_r602" "z0920_r602" "z0880_r602" ...
    %             "z0840_r602" "z0800_r602" "z0760_r602" "z0720_r602" ...
    %             "z0680_r602" "z0640_r602" "z0600_r602" "z0560_r602" ...
    %             "z0520_r602" "z0480_r602" "z0440_r602" "z0400_r602" ...
    %             "z0360_r602" "z0320_r602" "z0280_r602" "z0240_r602" ...
    %             "z0200_r602" "z0160_r602" "z0120_r602" "z0080_r602"];
% data_dirs = ["z0200_r602" "z0160_r602" "z0120_r602" "z0080_r602"];
close all;
data_dirs = ["z2033_r602"];
data_dirs = ["z0800_r602"];
outputDir = "OUTPUT/";

dir3c = "outputfiles_03c/";
env3c = load(dir3c + "equiv").param;

% グラフ
figure()
z = linspace(env3c.zmin, env3c.zmax, env3c.Nz);
r = linspace(env3c.rmin, env3c.rmax, env3c.Nr);
v = linspace(-20, 20, 21);

Nz = 2033;
Nr = 800;
zmin = -9.985000000000001e-01;
zmax = 9.985000000000001e-01;
rmin = 1.081500000000000e-01;
rmax = 0.8880;
delr = 9.747920133111480e-04;
delz = 9.827755905511811e-04;

zz = linspace(zmin, zmax, Nz);
rr = linspace(rmin, rmax, Nr);

for i = 1:length(data_dirs)
    vars = load(outputDir + data_dirs(i) + "/vars.mat");
    psi = vars.psi';

    % ０で埋める
    % A = size(vars.psi);
    % psi00 = zeros(env3c.Nr, env3c.Nz - A(1, 1));
    % psi = [vars.psi', psi00]';

    % 外側の磁束で埋める
    psiOutside = load("psiOutside.mat").A;
    padding = zeros(198, 2033);
    psi602 = [psi; padding];
    psi = psi602 + psiOutside';

    I = islocalmin(psi);
    I = find(I)
    [row, col] = ind2sub(size(psi), I);

    [M y] = min(min(psi));
    [M x] = min(psi(:, y));

    [a b] = islocalmin(min(psi));

    % プロット
    subplot(1, length(data_dirs), i)
    contour(rr, zz, psi' * 1000, 'r')
    hold on
    scatter(rr(row), zz(col));
    scatter(rr(x), zz(y));
    title(data_dirs(i))

end
