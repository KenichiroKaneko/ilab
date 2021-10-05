close all;

% 変数
select_phase = 6;
time_CCS = 9.45;
eddy_time = 0;
SENSio = 0;
LWALL = 0;
IUTST = 5;

IMAX = 512; % IMAX=513 (�d���f�[�^�폜)
MAXM = 10;
RSENS = 0.113150;
RWALL = 0.108150;
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

% 角度校正と読み込むデータ
foldername_TF = '180515';
shotnum_TF = '001';
foldername = '180515';
shotnum = '010';
% 角度校正用
% [time, Plasma_Current_TFshot, Coil_Current_TFshot, Flux_Loop_TFshot, Magnetic_Probe_TFshot, Magnetic_Probe_Lowpass_TFshot] = read_CCS_data(foldername_TF, shotnum_TF);
[time, Plasma_Current, Coil_Current, Flux_Loop, Magnetic_Probe, Magnetic_Probe_Lowpass] = read_CCS_data(foldername, shotnum);

% センサー位置の読み込みFLとMP
fid61 = fopen('CCS_FLXLP_sensor_position.txt', 'r');
fid68 = fopen('CCS_MP_sensor_position.txt', 'r');
textscan(fid61, '%s', 1, 'Delimiter', '\n');
textscan(fid68, '%s', 1, 'Delimiter', '\n');
% z_flxlp = zeros(1, 35);
% r_flxlp = ones(1, 35);
% z_mp = zeros(1, 40);
% r_mp = ones(1, 40);

% これなぞ
[Bz_EF_at_sensor_f, Psi_EF_at_sensor_f] = EF_calc_for_CCS_probe(r_flxlp, z_flxlp, 1, 19, EF_voltage);
[Bz_EF_at_sensor_b, Psi_EF_at_sensor_b] = EF_calc_for_CCS_probe(r_mp, z_mp, 1, 40, EF_voltage);

for i = 1:35
    temp = textscan(fid61, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    r_flxlp(1, i) = str2double(temp_m(1));
    z_flxlp(1, i) = str2double(temp_m(2));
end

% なにこれ
r_flxlp = r_flxlp .* 0.108785;

for i = 1:40
    temp = textscan(fid68, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    r_mp(1, i) = str2double(temp_m(1));
    z_mp(1, i) = str2double(temp_m(2));
end

TF_txt = zeros(40, 30000);
BZ_Distribution = zeros(1, MP_size(1));
BZ_normalization_value = zeros(1, MP_size(1));
EF_Distribution_MP = zeros(1, MP_size(1));

% 磁場センサー信号の整形、どこの時間で切り取るか
for i = 1:MP_size(1)
    TF_txt(i, :) = TF(1 + (i - 1) * 30000:i * 30000) * 0.001;
    BZ_t(i, :) = cumtrapz(time, Magnetic_Probe(i, :));
    BZ_t_Lowpass(i, :) = cumtrapz(Magnetic_Probe_Lowpass(i, :) - mean(Magnetic_Probe_Lowpass(i, 1:1000))) * 0.5e-6 * 242;
    % 角度校正用
    % BZ_t_TFshot(i, :) = cumtrapz(Magnetic_Probe_Lowpass_TFshot(i, :) - mean(Magnetic_Probe_Lowpass_TFshot(i, 1:1000))) * 0.5e-6 * 242;
    BZ_t_TFshot(i, :) = cumtrapz(Magnetic_Probe_Lowpass_TFshot(i, :)) * 0.5e-6 * 242;
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

% 磁束
% FLXLPの1~19chとそれ以外で違う操作をしている
% なぞい
% 多分インボードサイドが1~19ch
for i = 1:FL_size(1)

    if and(1 <= i, i <= 19)
        % 角度校正用
        % PSI_t(i, :) = (Flux_Loop(i, :) - Flux_Loop_TFshot(i, :)) / 202;
        PSI_t(i, :) = Flux_Loop(i, :) / 202;
        PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
        PSI_Distribution(1, i) = PSI_t(i, round(time_CCS / (0.5 * 0.001)));
        PSI_normalization_value(1, i) = var(PSI_t(i, 1:1000));
    else
        % 角度校正用
        % PSI_t(i, :) = (Flux_Loop(i, :) - Flux_Loop_TFshot(i, :)) / 37;
        PSI_t(i, :) = Flux_Loop(i, :) / 37;
        PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
        PSI_Distribution(1, i) = PSI_t(i, round(time_CCS / (0.5 * 0.001)));
        PSI_normalization_value(1, i) = var(PSI_t(i, 1:1000));
    end

end

PSI_Distribution = PSI_Distribution + Psi_EF_at_sensor_f;

figure()
plot(BZ_normalization_value)
figure()
plot(PSI_Distribution)

% reference
sensordata_B0 = fileread(['Sensor_B.txt']);
sensordata_B = strsplit(sensordata_B0, {'\n', '\t', '\r'});
sensornum_B = (length(sensordata_B) - 1) / 5 - 1;

for i = 1:sensornum_B
    BZ(i) = str2double(sensordata_B{9 + (i - 1) * 5});
    BR(i) = str2double(sensordata_B{10 + (i - 1) * 5});
end

sensordata_Flux0 = fileread(['Sensor_Flux.txt']);
sensordata_Flux = strsplit(sensordata_Flux0, {'\n', '\t', '\r'});
sensornum_Flux = (length(sensordata_Flux) - 1) / 5 - 1;

for i = 1:sensornum_Flux
    PSI(i) = str2double(sensordata_Flux{8 + (i - 1) * 5});
end

% KONDO
sensordata_Flux0 = fileread(['Parameters_FL_clockwise_180515010_t9450_EFkaizen.txt']);
sensordata_Flux = strsplit(sensordata_Flux0, {'\n', '\t', '\r'});
sensornum_Flux = (length(sensordata_Flux) - 1) / 5 - 1;

for i = 1:sensornum_Flux
    PSI_kondo(i) = str2double(sensordata_Flux{8 + (i - 1) * 5});
end

sensordata_B0 = fileread(['Parameters_MP_clockwise_180515010_t9450_EFkaizen.txt']);
sensordata_B = strsplit(sensordata_B0, {'\n', '\t', '\r'});
sensornum_B = (length(sensordata_B) - 1) / 5 - 1;

for i = 1:sensornum_B
    BZ_kondo(i) = str2double(sensordata_B{9 + (i - 1) * 5});
    BR_kondo(i) = str2double(sensordata_B{10 + (i - 1) * 5});
end

figure('Name', 'FLreference')
plot(PSI)
figure("Name", 'FLnama')
plot(PSI_Distribution(1, :))
figure('Name', 'FLkondo')
plot(PSI_kondo)

figure("Name", "MPreference")
plot(BZ)
figure("Name", "MPnama")
plot(BZ_Distribution(1, :))
figure("Name", "MPKondo")
plot(BZ_kondo)

fid61 = fopen('CCS_FLXLP_sensor_position.txt', 'r');
