% matファイルでの読み込み
vars = load("./CCS_dataset_gs/UTST_numel_2033" + "/merged.mat");
env = vars.env3c;

% プラズマ無しの容器内の磁束
psi602 = zeros(env.Nz, env.Nr_orig);
[env, psi602] = cal_psi(env, psi602);
psi602 = psi602.';

% プラズマ無しの容器外も含めた磁束
env.Nr = 800;
env.rmax = 0.8880;
psi800 = zeros(env.Nz, env.Nr);
[env, psi800] = cal_psi(env, psi800);
psi800 = psi800.';

% プラズマ有りの容器内の磁束
psi_orig = vars.psi;

% プラズマ有りの容器外も含めた磁束
padding = zeros(env.Nr - env.Nr_orig, env.Nz);
psi602 = [psi602; padding];
psi_orig = [psi_orig; padding];
whos
psi = psi_orig - psi602 + psi800;

br = vars.Br;
bz = vars.Bz;
z = linspace(env.zmin, env.zmax, env.Nz);
r = linspace(env.rmin, env.rmax, env.Nr);

% figure()
% % contour(psi)
% v = linspace(-40, 40, 121);
% contour(r, z, flipud(psi' * 1000), v, 'r')

% figure()
% contour(r, z, flipud(psi800' * 1000), v, 'r')
% hold on
% contour(r, z, flipud(psi602' * 1000), v, 'c')

% figure()
% contour(r, z, flipud((psi800 - psi602)' * 1000), v)

%% 実際のセンサー配置でデータを読み込む
sensorPosB = fileread("./CCS_temporary/CCS_MP_sensor_position.txt");
sensorPosB = strsplit(sensorPosB, {'\n', '\t', '\r'});
sensorPosFL = fileread("./CCS_temporary/CCS_FLXLP_sensor_position_3.txt");
sensorPosFL = strsplit(sensorPosFL, {'\n', '\t', '\r'});
SENSOR_TPRB.NUM = (length(sensorPosB) - 5) / 5;
SENSOR_FLXLP.NUM = (length(sensorPosFL) - 5) / 5;
% センサー配置をテキストデータから取得
SENSOR_TPRB.R = str2double(sensorPosB(6:5:5 * SENSOR_TPRB.NUM + 5));
SENSOR_TPRB.Z = str2double(sensorPosB(7:5:5 * SENSOR_TPRB.NUM + 5));
SENSOR_FLXLP.R = str2double(sensorPosFL(6:5:5 * SENSOR_FLXLP.NUM + 5));
SENSOR_FLXLP.Z = str2double(sensorPosFL(7:5:5 * SENSOR_FLXLP.NUM + 5));
% 磁場や磁束のデータをmatファイルから取得
for i = 1:SENSOR_TPRB.NUM
    [m, I] = min(abs(r - SENSOR_TPRB.R(i)));
    indexB_R(i) = I;
    [m, I] = min(abs(z -SENSOR_TPRB.Z(i)));
    indexB_Z(i) = I;
end

for i = 1:SENSOR_FLXLP.NUM
    [m, I] = min(abs(r - SENSOR_FLXLP.R(i)));
    indexFL_R(i) = I;
    [m, I] = min(abs(z -SENSOR_FLXLP.Z(i)));
    indexFL_Z(i) = I;
end

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

% for i = 1:SENSOR_FLXLP.NUM
%     % [m, I] = min(abs(R - SENSOR_FLXLP.R(i)) + abs(Z - SENSOR_FLXLP.Z(i)));
% end

% index_B = circshift(index_B, 31)
% index_FL = circshift(index_FL, 26)
% figure()
% hold on
% scatter(r(indexB_R), z(indexB_Z));
% % scatter(r(indexFL_R), z(indexFL_Z));
% scatter(SENSOR_TPRB.R, SENSOR_TPRB.Z, '*')

% for i = 1:SENSOR_TPRB.NUM
%     BR(i) = br(indexB_R(i), indexB_Z(i));
%     BZ(i) = bz(indexB_R(i), indexB_Z(i));
% end

% for i = 1:SENSOR_FLXLP.NUM
%     PSI(i) = psi(indexFL_R(i), indexFL_Z(i));
% end

% aaa = PSI(indexFL_R(:), indexFL_Z(:));
