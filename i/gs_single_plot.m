% data_dirs = ["outputfiles_03c/" "outputfiles_04c/" "outputfiles_05c/" "outputfiles_test/"];
% data_dirs = ["outputfiles_test/"];
% data_dirs = ["z2033_r602" ... 
%             "z1000_r602" "z0960_r602" "z0920_r602" "z0880_r602" ... 
%             "z0840_r602" "z0800_r602" "z0760_r602" "z0720_r602" ... 
%             "z0680_r602" "z0640_r602" "z0600_r602" "z0560_r602" ... 
%             "z0520_r602" "z0480_r602" "z0440_r602" "z0400_r602" ... 
%             "z0360_r602" "z0320_r602" "z0280_r602" "z0240_r602" ... 
%             "z0200_r602" "z0160_r602" "z0120_r602" "z0080_r602"];
data_dirs = ["z0200_r602" "z0160_r602" "z0120_r602" "z0080_r602"];
outputDir = "OUTPUT/";

dir3c = "outputfiles_03c/";
env3c = load(dir3c + "equiv").param;

% グラフ
z = linspace(env3c.zmin, env3c.zmax, env3c.Nz);
r = linspace(env3c.rmin, env3c.rmax, env3c.Nr);
v = linspace(-20, 20, 21);

for i = 1:length(data_dirs)

    % データの読み込み
    data_dir = outputDir + data_dirs(i) + "/";
    vars = load(data_dir + "vars");
    param = vars.param;

    % ０で埋める
    A = size(vars.psi);
    psi00 = zeros(env3c.Nr, env3c.Nz - A(1, 1));
    psi = [vars.psi', psi00]';

    % プロット
    subplot(1, length(data_dirs), i)
    contour(r, z, psi * 1000, v, 'r')
    title(data_dirs(i))

end
