% 片側の数値解を全体に拡張したものをプロットする
close all;
% data_dirs = ["z2033_r6022021" ...
    %             "z1000_r602", "z0960_r602", "z0920_r602", "z0880_r602" ...
    %             "z0840_r602", "z0800_r602", "z0760_r602", "z0720_r602" ...
    %             "z0680_r602", "z0640_r602", "z0600_r602", "z0560_r602" ...
    %             "z0520_r602", "z0480_r602", "z0440_r602", "z0400_r602" ...
    %             "z0360_r602", "z0320_r602", "z0280_r602", "z0240_r602" ...
    %             "z0200_r602", "z0160_r602", "z0120_r602", "z0080_r602"];
% data_dirs = ["z2033_r6022021", "z1000_r602", "z0960_r602", ];
% data_dirs = ["z0800_r602"];
% data_dirs = ["z2033_r6022021" "z1000_r602", "z0840_r602", ...
    %             "z0680_r602", "z0520_r602", "z0360_r602", "z0200_r602", ];
data_dirs = ["z1300_r6022021"];
gsOutputDir = "GSPLOT_OUTPUT/";

% 変数の定義
vars3c = load("OUTPUT/z2033_r602/vars");
env3c = vars3c.param;
zz = linspace(env3c.zmin, env3c.zmax, env3c.Nz);
rr = linspace(env3c.rmin, env3c.rmax, env3c.Nr);
figure()

for i = 1:length(data_dirs)
    vars = load(gsOutputDir + data_dirs(i) + "/merged.mat");
    psi = vars.psi;

    v = linspace(-20, 20, 11);
    subplot(1, length(data_dirs), i);
    contour(rr, zz, psi' * 1000, v, 'k', "LineWidth", 3);
    xlabel("R[m]");
    ylabel("Z[m]");

    % ylim([-0.4 0.4])
    hold on
    % [M1 y] = min(min(psi));
    % [M2 x] = min(psi(:, y));
    % scatter(rr(x), zz(y))

end

% プラズマ中心を配列に保存する
% 変数の定義
% vars3c = load("OUTPUT/z2033_r602/vars");
% env3c = vars3c.param;
% zz = linspace(env3c.zmin, 0, 1017);
% rr = linspace(env3c.rmin, env3c.rmax, env3c.Nr);

% for i = 1:length(data_dirs)
%     vars = load(gsOutputDir + data_dirs(i) + "/merged.mat");
%     psi = vars.psi;
%     psi(:, 1018:end) = [];
%     [M1 y] = min(min(psi));
%     [M2 x] = min(psi(:, y));
%     CCScenterR(i) = x;
%     CCScenterZ(i) = y;
%     % figure()
%     % contour(rr, zz, psi' * 1000, v, 'r');
%     % hold on
%     % [M1 y] = min(min(psi));
%     % [M2 x] = min(psi(:, y));
%     % scatter(rr(x), zz(y))
% end

% CCSrr = rr(CCScenterR);
% CCSzz = abs(zz(CCScenterZ));
% CCScenterRZ = [CCScenterR' CCScenterZ' CCSrr' CCSzz'];
% save("CCScenterRZ", "CCScenterRZ");
