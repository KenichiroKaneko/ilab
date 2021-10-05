close all;

% 変数
select_phase = 6;
time_CCS = 9.51;
eddy_time = 0;
SENSio = 0;
LWALL = 0;
IUTST = 5;

IMAX = 512; % IMAX=513 (�d���f�[�^�폜)
MAXM = 10;
RSENS = 0.113150;
RWALL = 0.108150;
% ***********************************************************************
%RS = zeros(1,nsemx);
%ZS = zeros(1,nsemx);
%ITYPE = zeros(1,nsemx);
%TET = zeros(1,nsemx);
% ***********************************************************************
%     =================================================================
% probe & flux-loop data(position data)
%     -----------------------------------------------------------------
% NTPB   : number of T-probes.
% NNPB   : number of N-probes.
% NAPB   : number of all(T & N) probes.
% NFLX   : number of Flux-Loops.
%     =================================================================
TPRB = zeros(1, 1024);
FLXLP = zeros(1, 1024);
%
% UTST �̕ǂ̊􉽍\��������
RSEC = [0.694D0, 0.694D0, 0.5985D0, 0.5985D0, 0.10815D0, 0.10815D0, ...
        0.10815D0, 0.5985D0, 0.5985D0, 0.694D0, 0.694D0];
ZSEC = [0.0D0, 0.285D0, 0.285D0, 0.9985D0, 0.9985D0, 0.0D0, ...
        -0.9985D0, -0.9985D0, -0.285D0, -0.285D0, 0.0D0];

foldername_TF = '180515';
shotnum_TF = '001';
[time, Plasma_Current_TFshot, Coil_Current_TFshot, Flux_Loop_TFshot, Magnetic_Probe_TFshot, Magnetic_Probe_Lowpass_TFshot] = read_CCS_data(foldername_TF, shotnum_TF);
% foldername = '180515';
% shotnum = '009';
foldername = '180515';
shotnum = '010';
[time, Plasma_Current, Coil_Current, Flux_Loop, Magnetic_Probe, Magnetic_Probe_Lowpass] = read_CCS_data(foldername, shotnum);

size(Plasma_Current)

Plasma_Current_minus_offset = Plasma_Current - mean(Plasma_Current(1000:2000));
plot(time, Plasma_Current)
hold on
% 100us�ňړ�����
figure('Name', 'Plasma Current smooth', 'NumberTitle', 'off')
Plasma_Current_minus_offset_smooth = smoothdata(Plasma_Current_minus_offset, 'gaussian', 200);
plot(time, Plasma_Current_minus_offset_smooth)
xlim([8.5, 10.5])

figure('Name', 'Plasma Current', 'NumberTitle', 'off')
subplot(5, 1, 5); plot(time, Plasma_Current_minus_offset, 'color', [0.7 0.7 0.7], 'DisplayName', 'Plasma current (raw)')
hold on
subplot(5, 1, 5); plot(time, Plasma_Current_minus_offset_smooth, 'r', 'DisplayName', 'Plasma current (smoothing)')
xlim([8.5, 10.5])
ylim([-100, 150])
xlabel('Time [ms]');
ylabel('Current [kA]');
legend({}, 'FontSize', 5, 'Location', 'eastoutside')

hold on
plot(time_CCS, Plasma_Current_minus_offset(time_CCS / (0.5 * 0.001)), 'o')

PF_Current_for_CCS = zeros(4, 1);
figure('Name', 'Coil Current', 'NumberTitle', 'off')
% smoothdata(Plasma_Current_minus_offset, 'gaussian', 200, 'DisplayName', 'Plasma Current smoothed')

for i = 1:4
    % �㉺�̃R�C���̓d���l�̕��ς�p����
    PF_Current_for_CCS(i, 1) = (Coil_Current(2 * i - 1, round(time_CCS / (0.5 * 0.001))) + Coil_Current(2 * i, round(time_CCS / (0.5 * 0.001)))) / 2;
    % �ۂ�
    PF_Current_for_CCS(i, 1) = round(PF_Current_for_CCS(i, 1), 5);
    subplot(4, 1, i); plot(time, Coil_Current(2 * i - 1, :), time_CCS, Coil_Current(2 * i - 1, time_CCS / (0.5 * 0.001)), 'o', 'DisplayName', 'PF%d (upper side)')
    subplot(5, 1, i); plot(time, Coil_Current(2 * i - 1, :), 'color', [0.7 0.7 0.7], 'DisplayName', 'PF#   current (upper side)')
    xlim([8.5, 10.5])
    ylim([-50, 100])
    hold on
    subplot(5, 1, i); plot(time, Coil_Current(2 * i, :), 'color', [0.9 0.9 0.9], 'DisplayName', 'PF#   current (lower side)')
    hold on
    subplot(5, 1, i); plot(time, smoothdata(Coil_Current(2 * i - 1, :), 'gaussian', 200), 'r', 'DisplayName', 'Smoothed PF#   current (upper side)')
    hold on
    subplot(5, 1, i); plot(time, smoothdata(Coil_Current(2 * i, :), 'gaussian', 200), 'b', 'DisplayName', 'Smoothed PF#   current (lower side)')
    hold on
    legend({}, 'FontSize', 5, 'Location', 'eastoutside')

end

fprintf('%d', PF_Current_for_CCS)

MP_size = size(Magnetic_Probe);
BZ_t = zeros(MP_size(1), MP_size(2));
BZ_t_TFshot = zeros(MP_size(1), MP_size(2));
BZ_t_Lowpass = zeros(MP_size(1), MP_size(2));
BZ_t_minusTF = zeros(MP_size(1), MP_size(2));

fileid = fopen('TF.txt', 'r');
TF = fscanf(fileid, '%f');

fileid3 = fopen('psi_EF_at_FL_clockwise.txt', 'r');
EF_FL = fscanf(fileid3, '%f');

% �Z���T�ʒu�ł�EF�̊�^���v�Z
fid61 = fopen('CCS_FLXLP_sensor_position.txt', 'r');
fid68 = fopen('CCS_MP_sensor_position.txt', 'r');
textscan(fid61, '%s', 1, 'Delimiter', '\n');
z_flxlp = zeros(1, 35);
r_flxlp = ones(1, 35);

for i = 1:35
    temp = textscan(fid61, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    r_flxlp(1, i) = str2double(temp_m(1));
    z_flxlp(1, i) = str2double(temp_m(2));
end

% �����t���b�N�X���[�v��r=0.108785m�ɐݒu
% なにこれ
r_flxlp = r_flxlp .* 0.108785;

z_mp = zeros(1, 40);
r_mp = ones(1, 40);
textscan(fid68, '%s', 1, 'Delimiter', '\n');

for i = 1:40
    temp = textscan(fid68, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    r_mp(1, i) = str2double(temp_m(1));
    z_mp(1, i) = str2double(temp_m(2));
end

[Bz_EF_at_sensor_f, Psi_EF_at_sensor_f] = EF_calc_for_CCS_probe(r_flxlp, z_flxlp, 1, 19, EF_voltage);
[Bz_EF_at_sensor_b, Psi_EF_at_sensor_b] = EF_calc_for_CCS_probe(r_mp, z_mp, 1, 40, EF_voltage);

TF_txt = zeros(40, 30000);
BZ_Distribution = zeros(1, MP_size(1));
BZ_normalization_value = zeros(1, MP_size(1));
EF_Distribution_MP = zeros(1, MP_size(1));

% 磁場センサー信号の整形
for i = 1:MP_size(1)
    TF_txt(i, :) = TF(1 + (i - 1) * 30000:i * 30000) * 0.001;
    BZ_t(i, :) = cumtrapz(time, Magnetic_Probe(i, :));
    BZ_t_Lowpass(i, :) = cumtrapz(Magnetic_Probe_Lowpass(i, :) - mean(Magnetic_Probe_Lowpass(i, 1:1000))) * 0.5e-6 * 242;
    BZ_t_TFshot(i, :) = cumtrapz(Magnetic_Probe_Lowpass_TFshot(i, :) - mean(Magnetic_Probe_Lowpass_TFshot(i, 1:1000))) * 0.5e-6 * 242;
    BZ_t_minusTF(i, :) = BZ_t_Lowpass(i, :) - BZ_t_TFshot(i, :);
    BZ_Distribution(1, i) = BZ_t_minusTF(i, round(time_CCS / (0.5 * 0.001)));
end

BZ_Distribution = BZ_Distribution + Bz_EF_at_sensor_b;
BZ_Distribution(30) = [];
BZ_normalization_value(30) = [];

FL_size = size(Flux_Loop);
PSI_t = zeros(FL_size(1), FL_size(2));
PSI_Distribution = zeros(1, FL_size(1));
PSI_normalization_value = zeros(1, FL_size(1));
EF_Distribution_FL = zeros(1, FL_size(1));

% 磁束センサー信号の整形
% FLXLPの1~19chとそれ以外で違う操作をしている
% なぞい
% 多分インボードサイドが1~19ch
for i = 1:FL_size(1)

    if and(1 <= i, i <= 19)
        PSI_t(i, :) = (Flux_Loop(i, :) - Flux_Loop_TFshot(i, :)) / 202;
        PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
        PSI_Distribution(1, i) = PSI_t(i, round(time_CCS / (0.5 * 0.001)));
        PSI_normalization_value(1, i) = var(PSI_t(i, 1:1000));
    else
        PSI_t(i, :) = (Flux_Loop(i, :) - Flux_Loop_TFshot(i, :)) / 37;
        PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
        PSI_Distribution(1, i) = PSI_t(i, round(time_CCS / (0.5 * 0.001)));
        PSI_normalization_value(1, i) = var(PSI_t(i, 1:1000));
    end

end

% (kondo)
figure('Name', 'fluxloop?')
subplot(1, 2, 1)
plot(PSI_normalization_value(1, 1:19).^(0.5))
hold on
plot([0, 20], [mean(PSI_normalization_value(1, 1:19).^(0.5)), mean(PSI_normalization_value(1, 1:19).^(0.5))])
subplot(1, 2, 2)
plot(PSI_normalization_value(1, 20:35).^(0.5))
hold on
plot([0, 20], [mean(PSI_normalization_value(1, 20:35).^(0.5)), mean(PSI_normalization_value(1, 20:35).^(0.5))])
mean(PSI_normalization_value(1, 1:19).^(0.5));
mean(PSI_normalization_value(1, 20:35).^(0.5));
% pause

% PSI_Distribution = PSI_Distribution + EF_Distribution_FL;
PSI_Distribution = PSI_Distribution + Psi_EF_at_sensor_f;

figure('Name', 'Poloidal Flux Distribution (all) including EF effect', 'NumberTitle', 'off')
PSI_Distribution(6) = [];
PSI_Distribution(35 - 1) = [];
plot(PSI_Distribution)

fid61 = fopen('Parameters_FL_clockwise_180515010_t8000_EFkaizen.txt', 'r');
fid68 = fopen('Parameters_MP_clockwise_180515010_t8000_EFkaizen.txt', 'r');
IMAX_b = 39;
IMAX_f = 33;

fid60 = fopen('@UTST_CheckWrite.txt', 'w');
%%%����fid62���v���n�̍��W�ɂ���
fid62 = fopen('@UTST_SenPos.txt', 'w');
fid63 = fopen('@UTST_VECTOR.txt', 'w');
fid64 = fopen('@UTST_WallGeom.txt', 'w');
fid65 = fopen('@UTST_CoilGeom.txt', 'w');
% ���C�f�[�^�ǂݍ���
textscan(fid61, '%s', 1, 'Delimiter', '\n'); % ��s�Ƃ΂�

% ***************************************************** %
%   Flux Loop                                           %
% ***************************************************** %
frewind(fid61);
textscan(fid61, '%s', 1, 'Delimiter', '\n'); %�������傤�Ƃ΂�
II = 0;
%
LGP = 1; % �����f�[�^���΂��������Ō���I�I�I�I17 15
LGPH = floor(LGP / 2) + 1;
LG = 0;

for I = 1:IMAX_f
    temp = textscan(fid61, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R = str2double(temp_m(1));
    Z = str2double(temp_m(2));
    %     PSI = str2double(temp_m(3));
    PSI = PSI_Distribution(1, I);

    LG = LG + 1;
    LLG = LG;

    if (LG == LGP) % 1�����Ń��Z�b�g����
        LG = 0;
    end

    if (LLG ~= LGPH) % 1������܂ł͓ǂ݂Ƃ΂�
        continue
    end

    if and(LWALL > 0, and (R < 0.12, abs(Z) < 0.9)); % �ǉ����̃Z���T�[����2016.9.5ushiki
        continue
    end

    E = abs((R - RSENS) / RSENS);
    II = II + 1;

    % センサー信号書き込み
    if (IUTST <= 9)
        RS(2 * II - 1) = R;
        RS(2 * II) = R;
        FLXLP(2 * II - 1) = PSI;
        FLXLP(2 * II) = PSI;
        %     end
        ZS(2 * II - 1) = Z;
        ZS(2 * II) = -Z;
        TET(2 * II - 1) = 0.0D0;
        TET(2 * II) = 0.0D0;
        ITYPE(2 * II - 1) = 0;
        ITYPE(2 * II) = 0;
        fprintf(fid60, '%d %d %d %d\r\n', 2 * II - 1, RS(2 * II - 1), ZS(2 * II - 1), FLXLP(2 * II - 1));
        fprintf(fid60, '%d %d %d %d\r\n', 2 * II, RS(2 * II), ZS(2 * II), FLXLP(2 * II));
    end

end

% �ߓ�
z_sp_max = 0.6;
z_sp_min = -0.6;
del_z_sp = 0.01;
zz = z_sp_min:del_z_sp:z_sp_max;
zz_pickup = [-0.2, -0.1, 0, 0.1, 0.2];
zz_idx = round((zz_pickup + 0.60) / 0.01 + 1);
ZS
% FLXLP_inboard_spline = spline(ZS(1:18), FLXLP(1:18), zz);
size(zz)
zz(zz_idx)

f_Bz_EF_z200 = fopen('180515009_9650_psi_r_z0.2.dat', 'r');
f_Bz_EF_z100 = fopen('180515009_9650_psi_r_z0.1.dat', 'r');
f_Bz_EF_z0 = fopen('180515009_9650_psi_r_z0.0.dat', 'r');
f_Bz_EF_z_nega_100 = fopen('180515009_9650_psi_r_z-0.1.dat', 'r');
f_Bz_EF_z_nega_200 = fopen('180515009_9650_psi_r_z-0.2.dat', 'r');

for i = 1:4
    temp = textscan(f_Bz_EF_z200, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1})
end

R_2d = zeros(1, 50)
psi = zeros(1, 50)

for i = 1:50
    temp = textscan(f_Bz_EF_z200, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1))
    psi(1, i) = str2double(temp_m(2))
end

plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(5));
Psi_z200 = psi + plus_inboard_flux;

for i = 1:4
    temp = textscan(f_Bz_EF_z100, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1})
end

R_2d = zeros(1, 50)
psi = zeros(1, 50)

for i = 1:50
    temp = textscan(f_Bz_EF_z100, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1))
    psi(1, i) = str2double(temp_m(2))
end

plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(4));
Psi_z100 = psi + plus_inboard_flux;

for i = 1:4
    temp = textscan(f_Bz_EF_z0, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
end

R_2d = zeros(1, 50);
psi = zeros(1, 50);

for i = 1:50
    temp = textscan(f_Bz_EF_z0, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1));
    psi(1, i) = str2double(temp_m(2));
end

plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(3));
Psi_z0 = psi + plus_inboard_flux;

for i = 1:4
    temp = textscan(f_Bz_EF_z_nega_100, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
end

R_2d = zeros(1, 50);
psi = zeros(1, 50);

for i = 1:50
    temp = textscan(f_Bz_EF_z_nega_100, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1));
    psi(1, i) = str2double(temp_m(2));
end

plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(2));
Psi_z_nega_100 = psi + plus_inboard_flux;

for i = 1:4
    temp = textscan(f_Bz_EF_z_nega_200, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
end

R_2d = zeros(1, 50);
psi = zeros(1, 50);

for i = 1:50
    temp = textscan(f_Bz_EF_z_nega_200, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1));
    psi(1, i) = str2double(temp_m(2));
end

plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(1));
Psi_z_nega_200 = psi + plus_inboard_flux;

if (IUTST <= 9)
    NFLX = 2 * II;
end

fprintf('%s %d\r\n', 'NFLX=', NFLX);
%
% ******************************************************
%   T-Probe
% ******************************************************
frewind(fid68); %
textscan(fid68, '%s', 1, 'Delimiter', '\n');
FACTR = 1.5D0;
II = 0;
%%C#####      LGB=1
LGB = 1; %21
LGBH = floor(LGB / 2) + 1;
LG = 0;
%%C####      MGB=1
MGB = 1; % �����ƂɃf�[�^���Ԉ����Ă��邩�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I
MGBH = floor(MGB / 2) + 1;
MG = 0;

for I = 1:IMAX_b
    temp = textscan(fid68, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R = str2double(temp_m(1));
    Z = str2double(temp_m(2));
    BR = str2double(temp_m(5)); %      READ(61,*),R,Z,PSI,BZ,BR
    BZ = BZ_Distribution(1, I);

    E = abs((R - RSENS) / RSENS);
    LG = LG + 1;
    LLG = LG;

    if (LG == LGB) % 1�����Ń��Z�b�g����
        LG = 0;
    end

    MG = MG + 1;
    MMG = MG;

    if (MG == MGB)
        MG = 0;
    end

    %
    %%C      IF(LLG.NE.LGBH) GOTO 200
    if and (E < 1.0E-5, MMG ~= MGBH)
        continue
    elseif and (E >= 1.0e-5, LLG ~= LGBH)
        continue

    elseif and(LWALL > 0, and (R < 0.12, abs(Z) < 0.9)); % �ǉ����̃Z���T�[����2016.9.5ushiki
        continue
    end

    II = II + 1;

    if (IUTST <= 9)
        RS(NFLX + 2 * II - 1) = R;
        RS(NFLX + 2 * II) = R;
        ZS(NFLX + 2 * II - 1) = Z;
        ZS(NFLX + 2 * II) = -Z;
        Z1 = Z; BR1 = BR;
        Z2 = -Z; BR2 = -BR;
        TET1 = atan2(BZ, BR1);
        TET2 = atan2(BZ, BR2);
        % ���͐��̐ڐ������̎�����E���I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I
        TET(NFLX + 2 * II - 1) = TET1; %  !
        TET(NFLX + 2 * II) = TET2; %  !�@
        ITYPE(NFLX + 2 * II - 1) = 1;
        ITYPE(NFLX + 2 * II) = 1;
        BB1 = BR1 * cos(TET1) + BZ * sin(TET1);
        BB2 = BR2 * cos(TET2) + BZ * sin(TET2);
        XBR1 = BB1 * cos(TET1);
        XBZ1 = BB1 * sin(TET1);
        XBR2 = BB2 * cos(TET2);
        XBZ2 = BB2 * sin(TET2);
        fprintf(fid60, '%d %d %d %d %d %d %d %d\r\n', R, Z1, BB1, TET1, XBR1, BR1, XBZ1, BZ);
        EPS1R = abs((XBR1 - BR1) / BR1);
        EPS1Z = abs((XBZ1 - BZ) / BZ);
        fprintf(fid60, '%d %d %d %d %d %d %d %d\r\n', R, Z2, BB2, TET2, XBR2, BR2, XBZ2, BZ);
        EPS2R = abs((XBR2 - BR2) / BR2);
        EPS2Z = abs((XBZ2 - BZ) / BZ);
        TPRB(2 * II - 1) = BB1;
        TPRB(2 * II) = BB2;
    end

end %200

if (IUTST <= 9)
    NTPB = 2 * II;
end

fprintf('%s %d\r\n', 'NTPB=', NTPB);
%
% N-Probe
%     Not defined.
NNPB = 0;
%%      NMAX=NFLX+NAPB
NMAX = NFLX + NTPB;
%

fid110 = fopen('SENPOS0.txt', 'w'); % 110
fid111 = fopen('SENPOS1.txt', 'w'); % 111
fid112 = fopen('SENPOS2.txt', 'w'); % 112

for I = 1:NMAX

    if (ITYPE(I) == 0)
        fprintf(fid110, '%d %d\r\n', RS(I), ZS(I));
    elseif (ITYPE(I) == 1)
        fprintf(fid111, '%d %d\r\n', RS(I), ZS(I));
    elseif (ITYPE(I) == 2)
        fprintf(fid112, '%d %d\r\n', RS(I), ZS(I));
    end

end

fclose(fid61);
fclose(fid60);
fclose(fid62);
fclose(fid63);
fclose(fid64);
fclose(fid65);
fclose(fid110);
fclose(fid111);
fclose(fid112);
%
%  *************************************************************************
%     Generation of CCS input data
%  *************************************************************************
%
fprintf('Generation of CCS input data ** START ***\n');
fid10 = fopen('CCSinput_UTST(temp).txt', 'w');
%
fprintf(fid10, '%s\r\n', '*');
fprintf(fid10, '%s\r\n', '*** CCS Test input for UTST generated in PreUTST ***');
fprintf(fid10, '%s\r\n', '** NTPB/NNPB/NFLX=(No. of T-/N-probes/Flux-Loops) **');
fprintf(fid10, '   %d     %d     %d\r\n', NTPB, NNPB, NFLX);
fprintf(fid10, '%s\r\n', '**  GETA (SSURF)');
fprintf(fid10, '  %s\r\n', '0.0000E+00');
fprintf(fid10, '%s\r\n', '* T-Probe');

for II = 1:NTPB
    fprintf(fid10, ' %d\r\n', TPRB(II));
end

fprintf(fid10, '%s\r\n', '* N-Probe');
fprintf(fid10, '%s\r\n', '* Flux-Loop');

for II = 1:NFLX
    fprintf(fid10, ' %d\r\n', FLXLP(II));
end

fprintf(fid10, '%s\r\n', '****** MINR * MAXR * MINZ * MAXZ ****');
fprintf(fid10, '  %s\r\n', '10   90  -100  100');
fprintf(fid10, '%s\r\n', '*********');
fprintf(fid10, '%s\r\n', '* ---�R�C���d���f�[�^�̕���---�P��[kA]');
fprintf(fid10, '%s\r\n', '* EF');
fprintf(fid10, '%s\r\n', '* PF#1');
fprintf(fid10, '%s\r\n', '* PF#2');
fprintf(fid10, '%s\r\n', '* PF#3');
fprintf(fid10, '%s\r\n', '* PF#4');

PF_Current_for_CCS = zeros(4, 1);

for i = 1:4
    % �㉺�̃R�C���̓d���l�̕��ς�p����
    PF_Current_for_CCS(i, 1) = (Coil_Current(2 * i - 1, round(time_CCS / (0.5 * 0.001))) + Coil_Current(2 * i, round(time_CCS / (0.5 * 0.001)))) / 2;
    %     subplot(4, 1, i); plot(time, Coil_Current(2 * i - 1, :), time_CCS, Coil_Current(2 * i - 1, time_CCS / (0.5 * 0.001)), 'o', 'DisplayName', 'PF%d (upper side)')
    %     hold on
    %     subplot(4, 1, i); plot(time, Coil_Current(2 * i, :), 'DisplayName', 'PF%d (lower side)')
    %     legend
end

fprintf('%d', PF_Current_for_CCS)

I_EF = -(0.849 * (1.19 * EF_voltage - 5.32) - 5.56);
% EF��200turn�R�C����6����ɂȂ��Ă���

fprintf(fid10, '%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10, '%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10, '%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10, '%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10, '%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10, '%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10, '%d\r\n', PF_Current_for_CCS(1, 1));
fprintf(fid10, '%d\r\n', PF_Current_for_CCS(2, 1));
fprintf(fid10, '%d\r\n', PF_Current_for_CCS(3, 1));
fprintf(fid10, '%d\r\n', PF_Current_for_CCS(4, 1));

fclose(fid10);
fprintf('Generation of CCS input data ** END ***\n');
