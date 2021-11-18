function [PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, ExtCOIL] = loadLizzieData(PARAM, CONFIG)
    % 生データをセンサーデータにするコード
    [BZ_Distribution, PSI_Distribution, ExtCOIL, PSI_normalization_value] = PreUTST(PARAM);
    sensorPosB = fileread("./CCS_temporary/CCS_MP_sensor_position_i.txt");
    sensorPosB = strsplit(sensorPosB, {'\n', '\t', '\r'});
    sensorPosFL = fileread("./CCS_temporary/CCS_FLXLP_sensor_position_i.txt");
    sensorPosFL = strsplit(sensorPosFL, {'\n', '\t', '\r'});
    SENSOR_TPRB.NUM = length(BZ_Distribution);
    SENSOR_FLXLP.NUM = length(PSI_Distribution);
    disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);
    disp(['Number of FLXLP = ' num2str(SENSOR_FLXLP.NUM)]);

    % SENSOR_TPRBを作成
    for i = 1:SENSOR_TPRB.NUM
        SENSOR_TPRB.R(i) = str2double(sensorPosB{6 + (i - 1) * 5});
        SENSOR_TPRB.Z(i) = str2double(sensorPosB{7 + (i - 1) * 5});
        SENSOR_TPRB.TET(i) = 1.570796326794897e+00;

        if SENSOR_TPRB.R(i) > 0.58
            SENSOR_TPRB.TET(i) = -1.570796326794897e+00;
        elseif SENSOR_TPRB.R(i) > 0.3
            SENSOR_TPRB.TET(i) = 0;
        end

        % SENSOR_TPRB.TET(i) = atan2(BZ_Distribution(i), 0);
        SENSOR_TPRB.TPRB(i) = abs(BZ_Distribution(i));
        SENSOR_TPRB.ITYPE(i) = 1; % 使っていない
        SENSOR_TPRB.TNPRB(i) = BZ_Distribution(i);
    end

    %% No NPRB
    SENSOR_NPRB.NUM = 0;
    SENSOR_NPRB.R = [];
    SENSOR_NPRB.Z = [];
    SENSOR_NPRB.TET = [];
    SENSOR_NPRB.NPRB = [];
    SENSOR_NPRB.ITYPE = [];

    % SENSOR_FLXLPを作成
    for i = 1:SENSOR_FLXLP.NUM
        SENSOR_FLXLP.R(i) = str2double(sensorPosFL{6 + (i - 1) * 5});
        SENSOR_FLXLP.Z(i) = str2double(sensorPosFL{7 + (i - 1) * 5});
        % SENSOR_FLXLP.FLXLP(i) = PSI_Distribution(i) / (pi * SENSOR_FLXLP.R(i)^2);
        SENSOR_FLXLP.FLXLP(i) = PSI_Distribution(i);
        % SENSOR_FLXLP.FLXLP(i) = PSI_normalization_value(i);
        % SENSOR_FLXLP.FLXLP(i) = 0;
        SENSOR_FLXLP.TET(i) = 1;
        SENSOR_FLXLP.ITYPE(i) = 0; % 使っていない
    end

    % CCSセンター位置の決定
    if CONFIG.DetermineCCSZPos
        RR = SENSOR_TPRB.R; ZZ = SENSOR_TPRB.Z;
        B = SENSOR_TPRB.TPRB;
        RR0index = RR < min(RR) + 0.1;
        B = B(RR0index); ZZ = ZZ(RR0index);
        ZZ = ZZ'; B = B';
        f = fit(ZZ, B, 'smoothingspline', 'SmoothingParam', 1);
        x = fnzeros(fnder(f.p));
        x = unique(x(:));
        y = fnval(f.p, x);
        matrix = [x y];
        matrix = sortrows(matrix, 2, 'descend');
        lmax = islocalmax(B);

        if CONFIG.ShowFig
            figure()
            hold on
            plot(ZZ, B);
            fnplt(f.p);
            plot(x, y, 'o');
            xlim([-1 1]);
            plot(ZZ(lmax), B(lmax), 'b*');
            view(90, 90);
        end

        if 1 %length(ZZ(lmax)) > 1
            PARAM.CCS = 2;
            kanekoZ = [0.14, -0.14];
            kanekoZ = [0.06, -0.06];

            for i = 1:2

                PARAM.R0(i) = 0.35;
                PARAM.Z0(i) = matrix(i, 1);
                % PARAM.Z0(i) = kanekoZ(i);
                PARAM.RR(i) = PARAM.RR(1);
                PARAM.CAPPER(i) = PARAM.CAPPER(1);
                PARAM.NE(i) = PARAM.NE(1);
                PARAM.MSEC(i, 1) = PARAM.MSEC(1, 1);
                PARAM.MSEC(i, 2) = PARAM.MSEC(1, 2);
                PARAM.MSEC(i, 3) = PARAM.MSEC(1, 3);
                PARAM.SOU(i) = PARAM.SOU(1);

            end

        end

        % CCS面の自動決定R方向
        if CONFIG.DetermineCCSRPos
            PARAM.R0 = CalcPlasmaCenter(PARAM, CONFIG, PARAM.Z0);
        end

    end

end

function [BZ_Distribution, PSI_Distribution, ExtCOIL, PSI_normalization_value] = PreUTST(PARAM)

    % function [BZ_normalization_value, PSI_normalization_value] = PreUTST(PARAM)
    %     =================================================================
    % probe & flux-loop data(position data)
    %     -----------------------------------------------------------------
    % NTPB   : number of T-probes.
    % NNPB   : number of N-probes.
    % NAPB   : number of all(T & N) probes.
    % NFLX   : number of Flux-Loops.
    %     =================================================================

    % Lizzieのデータから特定のショット・時間のFL、BZのデータを取得
    EF_voltage = 120;
    foldername = PARAM.foldername;
    shotnum = PARAM.shotnum;
    time_CCS = round(PARAM.time_CCS, 3);
    % round(time_CCS / (0.5 * 0.001))で30000点のデータの値になる

    %% Here, the geometory should be derived from input file
    % コイル電流の値は後で読み込む
    ExtCOIL.NUM = 10;
    ExtCOIL.NAME = ['EFL' 'EFU' 'PF1L' 'PF1U' 'PF2L' 'PF2U' 'PF3L' 'PF3U' 'PF4L' 'PF4U'];
    ExtCOIL.R = [0.80, 0.80, 0.20, 0.20, 0.665, 0.665, 0.750, 0.750, 0.685, 0.685];
    ExtCOIL.Z = [-1.07, 1.07, -1.10, 1.10, -0.80, 0.80, -0.675, 0.675, -0.50, 0.50];
    ExtCOIL.C = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    ExtCOIL.N = [200, 200, 8, 8, 3, 3, 8, 8, 3, 3];
    ExtCOIL.I = [0.28, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    foldername_TF = foldername;
    shotnum_TF = '001';
    [time, Plasma_Current_TFshot, Coil_Current_TFshot, Flux_Loop_TFshot, Magnetic_Probe_TFshot, Magnetic_Probe_Lowpass_TFshot] = read_CCS_data(foldername_TF, shotnum_TF, PARAM);
    [time, Plasma_Current, Coil_Current, Flux_Loop, Magnetic_Probe, Magnetic_Probe_Lowpass] = read_CCS_data(foldername, shotnum, PARAM);

    Plasma_Current_minus_offset = Plasma_Current - mean(Plasma_Current(1000:2000));
    Plasma_Current_minus_offset_smooth = smoothdata(Plasma_Current_minus_offset, 'gaussian', 200);

    % プラズマ電流・PF#1〜#4の上下のコイル電流のプロット
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
    plot(time_CCS, Plasma_Current_minus_offset(round(time_CCS / (0.5 * 0.001))), 'o')
    PF_Current_for_CCS = zeros(4, 1);

    for i = 1:4
        PF_Current_for_CCS(i, 1) = (Coil_Current(2 * i - 1, round(time_CCS / (0.5 * 0.001))) + Coil_Current(2 * i, round(time_CCS / (0.5 * 0.001)))) / 2;
        PF_Current_for_CCS(i, 1) = round(PF_Current_for_CCS(i, 1), 5);
        subplot(5, 1, i); plot(time, Coil_Current(2 * i - 1, :), time_CCS, Coil_Current(2 * i - 1, round(time_CCS / (0.5 * 0.001))), 'o', 'DisplayName', 'PF%d (upper side)')
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

    % コイル電流の読み込み
    for i = 1:8
        f = smoothdata(Coil_Current(i, :), 'gaussian', 200);

        % データ（NJ）はUpper→Lowerの順だが、CCSコードの内部ではLower→Upperの順で読み込んでいるので
        % インデックスに気をつけて読み込んでいる
        if rem(i, 2)
            j = i + 3;
        else
            j = i + 1;
        end

        ExtCOIL.I(j) = f(round(time_CCS / (0.5 * 0.001))) * 1;
    end

    ExtCOIL.I(1) =- (0.849 * (1.19 * EF_voltage - 5.32) - 5.56) / 1000;
    ExtCOIL.I(2) =- (0.849 * (1.19 * EF_voltage - 5.32) - 5.56) / 1000;

    MP_size = size(Magnetic_Probe);
    BZ_t = zeros(MP_size(1), MP_size(2));
    BZ_t_TFshot = zeros(MP_size(1), MP_size(2));
    BZ_t_Lowpass = zeros(MP_size(1), MP_size(2));
    BZ_t_minusTF = zeros(MP_size(1), MP_size(2));

    fid61 = fopen('./CCS_temporary/CCS_FLXLP_sensor_position_i.txt', 'r');
    fid68 = fopen('./CCS_temporary/CCS_MP_sensor_position_i.txt', 'r');
    textscan(fid61, '%s', 1, 'Delimiter', '\n');
    z_flxlp = zeros(1, 35);
    r_flxlp = ones(1, 35);

    for i = 1:35 % FluxLoopの位置を読み取っている
        temp = textscan(fid61, '%s', 1, 'Delimiter', '\n');
        temp_m = strsplit(temp{1}{1});
        r_flxlp(1, i) = str2double(temp_m(1));
        z_flxlp(1, i) = str2double(temp_m(2));
    end

    % ここ謎
    r_flxlp = r_flxlp .* 0.108785;

    z_mp = zeros(1, 40);
    r_mp = ones(1, 40);
    textscan(fid68, '%s', 1, 'Delimiter', '\n');

    for i = 1:40 % MPの位置を読み取っている
        temp = textscan(fid68, '%s', 1, 'Delimiter', '\n');
        temp_m = strsplit(temp{1}{1});
        r_mp(1, i) = str2double(temp_m(1));
        z_mp(1, i) = str2double(temp_m(2));
    end

    % 各センサー位置にEFコイルが作る磁場磁束を計算している？
    [Bz_EF_at_sensor_f, Psi_EF_at_sensor_f] = EF_calc_for_CCS_probe(r_flxlp, z_flxlp, 1, 35, EF_voltage);
    [Bz_EF_at_sensor_b, Psi_EF_at_sensor_b] = EF_calc_for_CCS_probe(r_mp, z_mp, 1, 40, EF_voltage);

    BZ_Distribution = zeros(1, MP_size(1));
    BZ_normalization_value = zeros(1, MP_size(1));
    EF_Distribution_MP = zeros(1, MP_size(1));

    NS = 242; % PUコイルの巻き数ｘ断面積

    for i = 1:MP_size(1)
        BZ_t(i, :) = cumtrapz(time, Magnetic_Probe(i, :));
        BZ_t_Lowpass(i, :) = cumtrapz(Magnetic_Probe_Lowpass(i, :) - mean(Magnetic_Probe_Lowpass(i, 1:1000))) * 0.5e-6 * NS;
        BZ_t_TFshot(i, :) = cumtrapz(Magnetic_Probe_Lowpass_TFshot(i, :) - mean(Magnetic_Probe_Lowpass_TFshot(i, 1:1000))) * 0.5e-6 * NS;
        BZ_t_minusTF(i, :) = BZ_t_Lowpass(i, :) - BZ_t_TFshot(i, :);
        BZ_Distribution(1, i) = BZ_t_minusTF(i, round(time_CCS / (0.5 * 0.001)));
    end

    BZ_Distribution = BZ_Distribution + Bz_EF_at_sensor_b;

    FL_size = size(Flux_Loop);
    PSI_t = zeros(FL_size(1), FL_size(2));
    PSI_Distribution = zeros(1, FL_size(1));
    PSI_normalization_value = zeros(1, FL_size(1));
    EF_Distribution_FL = zeros(1, FL_size(1));

    % 202と37はそれぞれ増幅器の係数
    FL_bai1 = 202;
    FL_bai2 = 37;
    FL_bai1 = 1;
    FL_bai2 = 1;

    figure('Name', 'FL sensor')
    hold on
    xline(round(time_CCS / (0.5 * 0.001)))

    for i = 1:FL_size(1)

        if and(1 <= i, i <= 19)
            PSI_t(i, :) = (Flux_Loop(i, :) - Flux_Loop_TFshot(i, :)) / FL_bai1;
            PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
            PSI_Distribution(1, i) = PSI_t(i, round(time_CCS / (0.5 * 0.001)));
            PSI_normalization_value(1, i) = var(PSI_t(i, 1:1000));
        else
            PSI_t(i, :) = (Flux_Loop(i, :) - Flux_Loop_TFshot(i, :)) / FL_bai2;
            PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
            PSI_Distribution(1, i) = PSI_t(i, round(time_CCS / (0.5 * 0.001)));
            PSI_normalization_value(1, i) = var(PSI_t(i, 1:1000));
        end

        plot(PSI_t(i, :))
    end

    legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35')

    % (kondo)
    figure('Name', 'PSI_normalization_value')
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

    figure('Name', 'PSI')
    hold on
    plot(PSI_Distribution)
    plot(Psi_EF_at_sensor_f)

    PSI_Distribution = PSI_Distribution + Psi_EF_at_sensor_f;
    % plot(PSI_Distribution2)

    % psi800 = load('psi800').psi800;
    % POS = load('POS').POS;

    % for i = 1:35
    %     Psi_EF_at_sensor_f(i) = psi800(POS.indexFL_Z(i), POS.indexFL_R(i));
    % end

    % plot(Psi_EF_at_sensor_f)
    % PSI_Distribution = PSI_Distribution .* PSI_normalization_value + Psi_EF_at_sensor_f;
    % plot(PSI_Distribution)
    % PSI_Distribution = PSI_Distribution;

    % figure('name', 'sensorData Bz')
    % plot(BZ_normalization_value)
    % hold on
    % plot(BZ_Distribution)
    % figure('name', 'sensorData flux')
    % plot(PSI_Distribution)
    % hold on
    % plot(PSI_normalization_value)
end

% need shot_no
function [time, Plasma_Current, Coil_Current, Flux_Loop, Magnetic_Probe, Magnetic_Probe_Lowpass] = read_CCS_data(foldername, shotnum, PARAM)

    lowpass_freq = 1000 * 1000;

    hdr_NJ = [foldername, '/NJ', shotnum, '.hdr'];
    hdr_NJ = "./utst/210906/NJ007.hdr"; % FLの増幅率をchごとに書いたheaderファイル
    dat_NJ = [foldername, '/NJ', shotnum, '.dat'];

    % ######################
    % ## READ HEADER DATA ##
    % ######################

    % read first line
    temp_NI = dlmread(hdr_NJ, '\t', [0 0 0 1]);
    ver_NI = temp_NI(1, 1);
    nch_NI = temp_NI(1, 2);
    % % READ headr
    fileID = fopen(hdr_NJ);
    header_NI = textscan(fileID, '%u %s %u %f64 %f64 %f64 %f64', 'HeaderLines', 1);
    fclose(fileID);
    chno = vertcat(header_NI{1});
    name = vertcat(header_NI{2});
    n = vertcat(header_NI{3});
    dt = vertcat(header_NI{4});
    range = vertcat(header_NI{6});
    gain = vertcat(header_NI{7});

    % ##########################
    % ###  READ BINARY DATA  ###
    % ##########################

    n = n(1);
    fid_NI = fopen(dat_NJ, 'r', 'b');
    data = horzcat(fread(fid_NI, [n, nch_NI], 'int16', 'b'));
    fclose(fid_NI);

    % #############################################
    % ###   CALIBLATE THE BINARY DATA           ###
    % ###  BY USING INFORMATION OF HEADER FILE  ###
    % #############################################

    time = zeros(1, 30000);
    dt = dt(1);
    I = 1:n;
    time(I) = double(I) * dt / double(1000);
    data_binary = zeros(n, nch_NI);

    for I = 1:nch_NI

        for J = 1:n
            data_binary(J, I) = range(I) / 32768 * gain(I) * data(J, I);
        end

    end

    % M�~N�z��@M�F�`�����l����ށAN�F�v���_�i15ms�A30000�_�j
    % ##### input plsma current current data #####
    Plasma_Current = zeros(1, 30000);
    Plasma_Current(1, :) = data_binary(:, 1);
    % figure('Name','Plasma Current','NumberTitle','off')
    % plot(time, Plasma_Current)

    % ##### input coil current data #####
    % Coil_Current = zeros(5, 30000);
    % figure('Name','Coil Current','NumberTitle','off')
    data_NJ = UTSTdataload('./', PARAM.date, PARAM.shotnum, 'NJ');

    for i = 1:length(data_NJ)
        % data_NJ(i).time(end - 24:end) = [];
        % data_NJ(i).raw(end - 24:end) = [];
        % data_NJ(i).sig(end - 24:end) = [];
        offset = mean(data_NJ(i).sig(50:round(450 / data_NJ(i).freq / 1e6)));
        data_NJ(i).sigint = cumtrapz(data_NJ(i).sig - offset) * data_NJ(i).freq;
    end

    for i = 2:9

        if i == 9 % 16�Ԗڂ�PF4L_feed���i�[����Ă���
            Coil_Current(8, :) = data_binary(:, 16);
            % Coil_Current(8, :) = data_NJ(16).sig;
            %         subplot(8, 1, 8); plot(time, Coil_Current(8, :))
            %         legend(name(16));
        else
            Coil_Current(i - 1, :) = data_binary(:, i);
            % Coil_Current(i - 1, :) = data_NJ(i).sig;
            %         subplot(8, 1, i - 1); plot(time, Coil_Current(i - 1, :))
            %         legend(name(i));
        end

    end

    % ##### input flux loop data #####
    Flux_Loop = zeros(35, 30000);
    % figure('Name','Flux Loop Signal (Gain multiplied (here, 1 or -1))','NumberTitle','off')

    index = [19, 20, 21, 22, 27, 28, 29, 31, 32, 33, 37, 38, 41, 42, 44, 45, 48, 50, 51, 53];

    for i = 18:53

        if (find(index == i))
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif i == 18 % 18�Ԗڂ�FL25���i�[����Ă���
            Flux_Loop(25, :) = data_binary(:, i);
            %         subplot(5, 7, 25); plot(time, Flux_Loop(25, :))
            %         legend(name(i));
        elseif i == 43
            continue
        else
            Flux_Loop(i - 18, :) = data_binary(:, i);
            %         subplot(5, 7, i - 18); plot(time, Flux_Loop(i - 18, :))
            %         legend(name(i));
        end

    end

    % % f1 = figure('Name', 'nama')
    % % f2 = figure('Name', 'FFT')
    % freq = fftshift((0:length(Flux_Loop(1, :)) - 1) * 2000000 / length(Flux_Loop(1, :))) - 1000000;

    % for i = 1:35
    %     % figure(f1)
    %     % hold on
    %     % plot(Flux_Loop(i, :))
    %     % hold off
    %     FL = fft(Flux_Loop(i, :));
    %     FL(freq > lowpass_freq) = 0;
    %     FL(freq < -lowpass_freq) = 0;
    %     Flux_Loop(i, :) = ifft(FL);
    %     % figure(f2)
    %     % hold on
    %     % plot(Flux_Loop(i, :))
    %     % hold off
    % end

    % ##### input magnetic probe data #####
    Magnetic_Probe = zeros(40, 30000);
    Magnetic_Probe_FT = zeros(40, 30000);
    Magnetic_Probe_FT_LP = zeros(40, 30000);
    Magnetic_Probe_Lowpass = zeros(40, 30000);
    freq = fftshift((0:length(Magnetic_Probe(1, :)) -1) * 2000000 / length(Magnetic_Probe(1, :))) - 1000000;
    % figure('Name','Magnetic Probe Signal (Gain multiplied (here, 1 or -1))','NumberTitle','off')

    index = [74, 79, 80, 81, 82, 84, 85, 86, 92];

    for i = 54:93

        if find(index == i)
            Magnetic_Probe(i - 53, :) = -data_binary(:, i);
        else
            Magnetic_Probe(i - 53, :) = data_binary(:, i);
        end

    end

    for i = 54:93
        %     Magnetic_Probe(i - 53, :) = data_binary(:, i);
        Magnetic_Probe_FT(i - 53, :) = fft(Magnetic_Probe(i - 53, :));
        Magnetic_Probe_FT_lp = Magnetic_Probe_FT(i - 53, :);
        Magnetic_Probe_FT_lp(freq > lowpass_freq) = 0;
        Magnetic_Probe_FT_lp(freq < -lowpass_freq) = 0;
        Magnetic_Probe_FT_LP(i - 53, :) = Magnetic_Probe_FT_lp;
        Magnetic_Probe_Lowpass(i - 53, :) = ifft(Magnetic_Probe_FT_LP(i - 53, :));
        %     subplot(5, 8, i - 53); plot(time, Magnetic_Probe_Lowpass(i - 53, :))
        %     legend(name(i));
    end

end

function [Bz, Psi] = EF_calc_for_CCS_probe(r, z, r_size, z_size, V)
    R_EF = 855 * 1e-3;
    z_EF = 1.05;
    Turn_EF = 200;
    mu0 = 4 * pi * 1e-7;
    V = 120;
    I =- (0.849 * (1.19 * V - 5.32) - 5.56);
    Psi = zeros(z_size, r_size);

    for i = 1:1
        alpha_u = sqrt((R_EF + r).^2 + (z - z_EF).^2);
        alpha_l = sqrt((R_EF + r).^2 + (z + z_EF).^2);
        k2_u = (4 * R_EF * r) ./ alpha_u.^2;
        k_u = sqrt(k2_u);
        k2_l = (4 * R_EF * r) ./ alpha_l.^2;
        k_l = sqrt(k2_l);
        m_u = k2_u;
        m_l = k2_l;
        [K_u, E_u] = ellipke(m_u);
        [K_l, E_l] = ellipke(m_l);
        Bz_upper = ((mu0 * Turn_EF * I) ./ (2 * pi)) .* (1 ./ alpha_u) .* (K_u + ((R_EF^2 - r.^2 - (z - z_EF).^2) ./ ((R_EF - r).^2 + (z - z_EF).^2)) .* E_u);
        Bz_lower = ((mu0 * Turn_EF * I) ./ (2 * pi)) .* (1 ./ alpha_l) .* (K_l + ((R_EF^2 - r.^2 - (z + z_EF).^2) ./ ((R_EF - r).^2 + (z + z_EF).^2)) .* E_l);
        Bz = Bz_upper + Bz_lower;

        A_phai_u = (mu0 * Turn_EF * I) ./ (pi) .* sqrt(R_EF ./ r) .* (1 ./ k_u) .* ((1 - k_u.^2 ./ 2) .* K_u - E_u);
        A_phai_l = (mu0 * Turn_EF * I) ./ (pi) .* sqrt(R_EF ./ r) .* (1 ./ k_l) .* ((1 - k_l.^2 ./ 2) .* K_l - E_l);

        psi_u = 2 * pi * r .* A_phai_u;
        psi_l = 2 * pi * r .* A_phai_l;
        Psi = psi_u + psi_l;

    end

end
