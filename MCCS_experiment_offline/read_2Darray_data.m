% need shot_no
function [Psi_2D, Bz_distribution_plus_EF, Psi_2D_sp] = read_2Darray_data(foldername, shotnum, t, EF_voltage)

    if 0
        %%%%%%%%%% test %%%%%%%%%%
        close all
        foldername = '190812';
        shotnum = '019';
        t = 9.45
        EF_voltage = 120
        %%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    hdr_NI = [foldername + '/NI' + shotnum + '.hdr'];
    dat_NI = [foldername + '/NI' + shotnum + '.dat'];
    hdr_ICS = [foldername + '/ICS' + shotnum + '.hdr'];
    dat_ICS = [foldername + '/ICS' + shotnum + '.dat'];
    angle = [foldername + '/angle' + '.txt'];
    fileID = fopen(angle);
    Angle_file = textscan(fileID, '%u %s %u %f64 %f64 %f64 %f64', 'HeaderLines', 0);
    fclose(fileID);
    Angle = vertcat(Angle_file{2});
    Angle = str2double(Angle);
    Angle_rad = Angle * pi / 180;
    figure('Name', 'Angle of coil tilting', 'NumberTitle', 'off')
    plot(Angle_rad)
    %ver_ICS = temp_ICS(1,1);
    %nch_ICS = temp_ICS(1,2);

    % ######################
    % ## READ HEADER DATA ##
    % ######################

    % read first line
    temp_NI = dlmread(hdr_NI, '\t', [0 0 0 1]);
    ver_NI = temp_NI(1, 1);
    nch_NI = temp_NI(1, 2);
    temp_ICS = dlmread(hdr_ICS, '\t', [0 0 0 1]);
    ver_ICS = temp_ICS(1, 1);
    nch_ICS = temp_ICS(1, 2);
    % % READ header
    fileID = fopen(hdr_NI);
    header_NI = textscan(fileID, '%u %s %u %f64 %f64 %f64 %f64', 'HeaderLines', 1);
    fclose(fileID);
    fileID = fopen(hdr_ICS);
    header_ICS = textscan(fileID, '%u %s %u %f64 %f64 %f64 %f64', 'HeaderLines', 1);
    fclose(fileID);

    chno = vertcat(header_NI{1});
    name_NI = vertcat(header_NI{2});
    n_NI = vertcat(header_NI{3});
    dt_NI = vertcat(header_NI{4});
    range_NI = vertcat(header_NI{6});
    gain_NI = vertcat(header_NI{7});
    % ##########################
    % ###  READ BINARY DATA (NI)  ###
    % ##########################
    n_NI = n_NI(1);
    fid_NI = fopen(dat_NI, 'r', 'b');
    data_NI = horzcat(fread(fid_NI, [n_NI, nch_NI], 'int16', 'b'));
    fclose(fid_NI);
    % #############################################
    % ###   CALIBLATE THE BINARY DATA           ###
    % ###  BY USING INFORMATION OF HEADER FILE  ###
    % #############################################
    time = zeros(1, 30000);
    dt_NI = dt_NI(1);
    I = 1:n_NI;
    time(I) = double(I) * dt_NI / double(1000);
    data_binary_NI = zeros(n_NI, nch_NI);

    for I = 1:nch_NI;

        for J = 1:n_NI;
            data_binary_NI(J, I) = range_NI(I) / 32768 * gain_NI(I) * data_NI(J, I);
        end

    end

    chno = vertcat(header_ICS{1});
    name_ICS = vertcat(header_ICS{2});
    n_ICS = vertcat(header_ICS{3});
    dt_ICS = vertcat(header_ICS{4});
    range_ICS = vertcat(header_ICS{6});
    gain_ICS = vertcat(header_ICS{7});
    % ##########################
    % ###  READ BINARY DATA (ICS)  ###
    % ##########################
    n_ICS = n_ICS(1);
    fid_ICS = fopen(dat_ICS, 'r', 'b');
    data_ICS = horzcat(fread(fid_ICS, [n_ICS, nch_ICS], 'int16', 'b'));
    fclose(fid_ICS);
    % #############################################
    % ###   CALIBLATE THE BINARY DATA           ###
    % ###  BY USING INFORMATION OF HEADER FILE  ###
    % #############################################
    time = zeros(1, 30000);
    dt_ICS = dt_ICS(1);
    I = 1:n_ICS;
    time(I) = double(I) * dt_ICS / double(1000);
    data_binary_ICS = zeros(n_ICS, nch_ICS);

    for I = 1:nch_ICS

        for J = 1:n_ICS
            data_binary_ICS(J, I) = range_ICS(I) / 32768 * gain_ICS(I) * data_ICS(J, I);
        end

    end

    disp(name_ICS);
    disp(gain_ICS);

    % ICS��2~82��Bt�v��
    % NI��1~81��Bz�v��
    % for i = 1:8
    %     subplot(2, 4, i); plot(time, data_binary_ICS(1 + 30000 * (i -1) : 30000 * i))
    %     legend(name_ICS(i))
    % %     xlim([8.5, 10.5])
    % %     ylim([-5, 5])
    % %     if i == 54
    % %         xlabel('Time [us]');
    % %         ylabel('Axial Magnetic Field [T]');
    % end

    %����f�[�^�Z�b�g���쐬
    Bz = zeros(81, 30000);
    Bz_1 = zeros(81, 30000);
    Bt = zeros(81, 30000);

    % size(data_binary_NI)
    % size(data_binary_ICS)
    %

    for i = 1:81
        rawsig_NI = data_binary_NI(:, i);
        rawsig_ICS = data_binary_ICS(:, i + 1);
        Bz_1(i, :) = cumtrapz(rawsig_NI) * 0.5e-6;
        Bz(i, :) = cumtrapz(rawsig_NI - mean(rawsig_NI(1:1000))) * 0.5e-6;
        Bt(i, :) = cumtrapz(rawsig_ICS - mean(rawsig_ICS(1:1000))) * 0.5e-6;
    end

    % figure('Name','Measured Troidal Field','NumberTitle','off')
    % for j = 1 : 9
    %     for i = 1 : 9
    %         subplot(9, 9, 9 * (j - 1) + i); plot(time, Bt(9 * (j - 1) + i, :))
    % %         legend(name_ICS(i))
    %     end
    % end
    %
    % figure('Name','Measured Axial Field','NumberTitle','off')
    % for j = 1 : 9
    %     for i = 1 : 9
    %         subplot(9, 9, 9 * (j - 1) + i); plot(time, Bz(9 * (j - 1) + i, :))
    %         hold on
    %         subplot(9, 9, 9 * (j - 1) + i); plot(time, Bz_1(9 * (j - 1) + i, :))
    % %         legend(name_ICS(i))
    %     end
    % end
    %
    %

    % t = 9.51;
    % t = 9.53;
    % t = 9.65;

    for i = 1:81
        Bt_real(i, :) = Bz(i, :) * sin(Angle_rad(i)) + Bt(i, :) * cos(Angle_rad(i));
        Bz_real(i, :) = Bz(i, :) * cos(Angle_rad(i)) - Bt(i, :) * sin(Angle_rad(i));
    end

    Bz_distribution = zeros(9, 9);
    figure('Name', 'Real Axial Field', 'NumberTitle', 'off')

    for j = 1:9

        for i = 1:9
            %         subplot(9, 9, 9 * (j - 1) + i); plot(time, Bz_real(9 * (j - 1) + i, :))
            %         hold on
            BZ_real = Bz_real(9 * (j - 1) + i, :);
            Bz_distribution(j, i) = BZ_real(int16(t / (0.5 * 0.001)));
            %         subplot(9, 9, 9 * (j - 1) + i); plot(t, BZ_real(t / (0.5 * 0.001)), 'o')
            %         legend(name_ICS(i))
        end

    end

    r = [0.11, 0.18, 0.25, 0.32, 0.39, 0.46, 0.53, 0.60, 0.67];
    z = [0.23, 0.18, 0.12, 0.06, 0.00, -0.06, -0.12, -0.18, -0.23];
    % 2D�A���C�ɏ��ׂ�EF��^���Z�o
    [Bz_EF, Psi_EF] = EF_calc(r, z, 9, 9, EF_voltage);

    Bz_distribution_plus_EF = zeros(9, 9);
    R = zeros(9, 9);

    for i = 1:9
        R(i, :) = r;
    end

    % figure('Name','Bz radial distribution','NumberTitle','off')
    for i = 1:9
        Bz_distribution_plus_EF(i, :) = Bz_distribution(i, :) + Bz_EF(i, :);
        %     plot(r, Bz_distribution(i, :))
        %     hold on
        %     plot(r, Bz_distribution_plus_EF(i, :))
        %     hold on
    end

    % plot(z, Bz_distribution_plus_EF(:, 1))
    % hold on
    %

    zz = -0.23:0.01:0.23;
    zz = flip(zz);
    %
    % Bz_distribution_plus_EF_sp(1, :) = spline(z(2 : 8), Bz_distribution_plus_EF(2 : 8), zz);
    % plot(zz, Bz_distribution_plus_EF_sp)
    % plot(zz(1), Bz_distribution_plus_EF_sp(1), 'o')
    % plot(zz(47), Bz_distribution_plus_EF_sp(47), 'o')
    %
    % zz_ind = (zz_pick - (-0.23)) / 0.01 + 1;

    % �s���`�����l������
    % Bz(z)���z�ŃX�v���C��
    for i = 1:9
        Bz_distribution_plus_EF_sub = Bz_distribution_plus_EF(:, i);

        if i == 1
            z(1) = []
            z(9 - 1) = []
            Bz_distribution_plus_EF_sub(1) = []
            Bz_distribution_plus_EF_sub(9 - 1) = []
            Bz_distribution_plus_EF_sp = spline(z, Bz_distribution_plus_EF_sub, zz);
            zz_pick = 0.23;
            zz_idx = int8((zz_pick - 0.23) / (-0.01) + 1);
            Bz_distribution_plus_EF(1, i) = Bz_distribution_plus_EF_sp(zz_idx);
            zz_pick = -0.23;
            zz_idx = int8((zz_pick - 0.23) / (-0.01) + 1);
            Bz_distribution_plus_EF(9, i) = Bz_distribution_plus_EF_sp(zz_idx);
            z = [0.23, 0.18, 0.12, 0.06, 0.00, -0.06, -0.12, -0.18, -0.23];
            plot(z, Bz_distribution_plus_EF(:, i))
            hold on
        elseif i == 2
            z(8) = []
            Bz_distribution_plus_EF_sub(8) = []
            Bz_distribution_plus_EF_sp = spline(z, Bz_distribution_plus_EF_sub, zz);
            zz_pick = -0.18;
            zz_idx = int8((zz_pick - 0.23) / (-0.01) + 1);
            Bz_distribution_plus_EF(8, i) = Bz_distribution_plus_EF_sp(zz_idx);
            z = [0.23, 0.18, 0.12, 0.06, 0.00, -0.06, -0.12, -0.18, -0.23];
            plot(z, Bz_distribution_plus_EF(:, i))
            hold on
        elseif i == 3
            plot(z, Bz_distribution_plus_EF(:, i))
            hold on
        elseif i == 4
            plot(z, Bz_distribution_plus_EF(:, i))
            hold on
        elseif i == 5
            plot(z, Bz_distribution_plus_EF(:, i))
            hold on
        elseif i == 6
            plot(z, Bz_distribution_plus_EF(:, 6))
            hold on
        elseif i == 7
            z(1) = []
            Bz_distribution_plus_EF_sub(1) = []
            Bz_distribution_plus_EF_sp(1, :) = spline(z, Bz_distribution_plus_EF_sub, zz);
            zz_pick = 0.23;
            zz_idx = int8((zz_pick - 0.23) / (-0.01) + 1);
            Bz_distribution_plus_EF(1, i) = Bz_distribution_plus_EF_sp(zz_idx);
            z = [0.23, 0.18, 0.12, 0.06, 0.00, -0.06, -0.12, -0.18, -0.23];
            plot(z, Bz_distribution_plus_EF(:, i))
            hold on
        elseif i == 8
            plot(z, Bz_distribution_plus_EF(:, 8))
            hold on
        elseif i == 9
            plot(z, Bz_distribution_plus_EF(:, 9))
        end

    end

    %
    % %�s���`�����l�����X�v���C����Ԃœ��}����
    % rr = 0.11 : 0.01 : 0.67;
    % rr_size = size(rr);
    % Bz_distribution_plus_EF_sp = zeros(9, rr_size(2));
    % r(1) = [];
    % Bz_distribution_plus_EF_sub = Bz_distribution_plus_EF(1, :)
    % Bz_distribution_plus_EF_sub(1) = []
    % r(7 - 1) = [];
    % Bz_distribution_plus_EF_sub(7 - 1) = []
    % r(9 - 2) = [];
    % Bz_distribution_plus_EF_sub(9 - 2) = []

    % �Z�o����Bz2�����v���t�@�C����9�~8(0.11 <= r <= 0.60)
    Bz_2D = zeros(9, 8);
    figure('Name', 'Check Bz', 'NumberTitle', 'off')

    for i = 1:9
        Bz_distribution_plus_EF_sub2 = Bz_distribution_plus_EF(i, :);
        Bz_2D(i, :) = Bz_distribution_plus_EF_sub2(1:8);
        plot(r(1:8), Bz_2D(i, :))
        hold on
        %
    end

    %

    [time, Plasma_Current, Coil_Current, Flux_Loop, Magnetic_Probe, Magnetic_Probe_Lowpass] = read_CCS_data(foldername, shotnum);
    FL_size = size(Flux_Loop);
    PSI_t = zeros(FL_size(1), FL_size(2));
    PSI_Distribution = zeros(1, FL_size(1));
    EF_Distribution_FL = zeros(1, FL_size(1));
    % FL signal [V] = Gain * 1 / (RC) * (integral((d��/dt)dt))
    % Gain : 100 (inboard), 10 (outboard)
    % RC : 510 (inboard), 470 (outboard)

    fid61 = fopen('Parameters_FL_clockwise_180515010_t8000_EFkaizen.txt', 'r');
    textscan(fid61, '%s', 1, 'Delimiter', '\n');
    z_flxlp = zeros(1, 18);
    r_flxlp = ones(1, 18);

    for i = 1:18
        temp = textscan(fid61, '%s', 1, 'Delimiter', '\n');
        temp_m = strsplit(temp{1}{1});
        z_flxlp(1, i) = str2double(temp_m(2));
    end

    % �����t���b�N�X���[�v��r=0.108785m�ɐݒu
    r_flxlp = r_flxlp .* 0.108785;

    for i = 1:19
        %     PSI_t(i, :) = (510. * 0.001 / 100.) * Flux_Loop(i, :);
        PSI_t(i, :) = Flux_Loop(i, :) / 202;
        PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
        subplot(5, 7, i); plot(time, PSI_t(i, :), 'DisplayName', 'inboard')
        legend
        PSI_Distribution(1, i) = PSI_t(i, int16(t / (0.5 * 0.001)));
    end

    [Bz_EF, Psi_EF] = EF_calc_for_CCS_probe(r_flxlp, z_flxlp, 1, 18, EF_voltage);

    figure('Name', 'Poloidal Flux Distribution (inboard)', 'NumberTitle', 'off')
    PSI_Distribution(6) = [];
    plot(z_flxlp, PSI_Distribution(1:18))
    hold on
    plot(z_flxlp, Psi_EF)
    hold on
    plot(z_flxlp, PSI_Distribution(1:18) + Psi_EF)

    PSI_Distribution_plus_EF = PSI_Distribution(1:18) + Psi_EF;

    z_sp_max = 0.6;
    z_sp_min = -0.6;
    del_z_sp = 0.01;
    zz = z_sp_min:del_z_sp:z_sp_max; % zz = 0.01 * (n - 1) + (-0.60)
    % 2�����A���C�f�[�^�ɑ������킹�鎥���𒊏o�@���@FLXLP_inboard_spline
    % zz_pickup = 0.0;
    zz_pickup = [-0.23, -0.18, -0.12, -0.06, -0.0, 0.06, 0.12, 0.18, 0.23];
    zz_idx = round((zz_pickup + 0.60) / 0.01 + 1);
    FLXLP_inboard_spline = spline(z_flxlp, PSI_Distribution_plus_EF, zz);

    figure('Name', 'Poloidal Flux', 'NumberTitle', 'off')
    plot(z_flxlp, PSI_Distribution_plus_EF)
    hold on
    plot(zz, FLXLP_inboard_spline, 'r-')
    hold on
    plot(zz(zz_idx), FLXLP_inboard_spline(zz_idx), 'ro')

    initial_flux_for_2Darray = FLXLP_inboard_spline(zz_idx);

    % �����̌v�Z
    Psi_2D = zeros(9, 8);
    % �C���{�[�h���t���b�N�X���[�v�̌v���l�������l�ɗp����
    Psi_2D_initial = ones(9, 8);

    for i = 1:9
        Psi_2D_initial(i, :) = Psi_2D_initial(i, :) .* initial_flux_for_2Darray(i);
        %     Psi_2D_initial(i, :) = Psi_2D_initial(i, :) .* 0;
    end

    figure('Name', 'Psi radial distribution', 'NumberTitle', 'off')

    for i = 1:9
        ker = r(1:8) .* Bz_2D(i, :);
        Psi_2D(i, :) = Psi_2D_initial(i, :) + 2 * pi * cumtrapz(r(1:8), ker);
        %     plot(r(1 : 8), Psi_2D(i, :))
        plot(r(1:8), Bz_2D(i, :))
        hold on
        plot(r(1:8), ker)
        hold on
        %
    end

    % Psi_2D���X�v���C����3���⊮����i�������j�������i190304�j
    size(Psi_2D)
    % R�����ɕ��
    r = [0.11, 0.18, 0.25, 0.32, 0.39, 0.46, 0.53, 0.60];
    rr = 0.11:0.01:0.60;
    Psi_2D_sp = zeros(47, 50);

    for i = 1:9
        Psi_2D(i, :)

        if i == 1
            Psi_2D_sp(1, :) = spline(r, Psi_2D(i, :), rr);
            Psi_2D_sp_sub(i, :) = spline(r, Psi_2D(i, :), rr);
        elseif i == 9
            Psi_2D_sp(42 + 5, :) = spline(r, Psi_2D(i, :), rr);
            Psi_2D_sp_sub(i, :) = spline(r, Psi_2D(i, :), rr);
        else
            Psi_2D_sp(6 + 6 * (i - 2), :) = spline(r, Psi_2D(i, :), rr);
            Psi_2D_sp_sub(i, :) = spline(r, Psi_2D(i, :), rr);
        end

        %     figure('Name','Psi dist. in-board','NumberTitle','off')
        %     plot(r, Psi_2D(i, :), 'o')
        %     hold on
        %     plot(rr, Psi_2D_sp(i, :))
        %
    end

    size(Psi_2D_sp_sub)
    size_Psi_2D_sp = size(Psi_2D_sp);
    z = [0.23, 0.18, 0.12, 0.06, 0.00, -0.06, -0.12, -0.18, -0.23];
    zz = -0.23:0.01:0.23
    zz = sort(zz, 'descend')

    for i = 1:size_Psi_2D_sp(2)
        Psi_2D_sp(:, i) = spline(z, Psi_2D_sp_sub(:, i), zz);
    end

    min_Psi_2D_sp = min(min(Psi_2D_sp))
    max_Psi_2D_sp = max(max(Psi_2D_sp))
    del_contour_sp = (max_Psi_2D_sp - min_Psi_2D_sp) / 50
    figure('Name', 'Poloidal Flux (spline ver.)', 'NumberTitle', 'off')
    contourf(rr, zz, Psi_2D_sp, 30, 'ShowText', 'on')
    hold on
    contour(rr, zz, Psi_2D_sp, max(Psi_2D_sp(:, 1)):del_contour_sp * 0.001:max(Psi_2D_sp(:, 1)) + del_contour_sp * 0.001, 'r-', 'LineWidth', 2)

    figure('Name', 'Psi dist. in-board', 'NumberTitle', 'off')
    plot(Psi_2D(:, 1))
    max(Psi_2D(:, 1))
    % contourf(r(1 : 8), z, Psi_2D, 20)
    min_Psi_2D = min(min(Psi_2D))
    max_Psi_2D = max(max(Psi_2D))
    del_contour = (max_Psi_2D - min_Psi_2D) / 50
    %
    figure('Name', 'Psi 2D', 'NumberTitle', 'off')
    % contourf(r(1 : 8), z, Psi_2D, 15, 'ShowText', 'on')
    contourf(r(1:8), z, Psi_2D, min_Psi_2D:del_contour:max_Psi_2D, 'ShowText', 'on')
    hold on
    contour(r(1:8), z, Psi_2D, max(Psi_2D(:, 1)):del_contour * 0.001:max(Psi_2D(:, 1)) + del_contour * 0.001, 'r-', 'LineWidth', 2)
    c = colorbar;
    % max(Psi_2D(: , 1)) : del_contour * 0.001 : max(Psi_2D(: , 1))
    % caxis([-0.000, 0.007]);
    % caxis([-0.001 0.003]);
    xlabel({'r [m]'});
    ylabel({'z [m]'});
    axis equal
    xlim([0.1, 0.6])
    hold on

    for i = 1:9
        z_2d_geom = ones(1, 8) * z(i);
        plot(r(1:8), z_2d_geom, 'kx')
        hold on
    end

    plot(r_flxlp(7:11), z_flxlp(7:11), 'ko')
    hold on
    plot([0.1085, 0.1085], [-0.23, 0.23], 'r--')
    hold on
    fid61 = fopen('Parameters_MP_clockwise_180515010_t8000_EFkaizen.txt', 'r');
    textscan(fid61, '%s', 1, 'Delimiter', '\n');
    z_mp = zeros(1, 18);
    r_mp = zeros(1, 18);

    for i = 1:18
        temp = textscan(fid61, '%s', 1, 'Delimiter', '\n');
        temp_m = strsplit(temp{1}{1});
        r_mp(1, i) = str2double(temp_m(1));
        z_mp(1, i) = str2double(temp_m(2));
    end

    plot(r_mp(8:11), z_mp(8:11), 'k^')

    % Br���Z�o
    for j = 1:size_Psi_2D_sp(1)

        for i = 1:size_Psi_2D_sp(2) - 1
            dpsi_dr(i, :) = (Psi_2D_sp(j, i + 1) - Psi_2D_sp(j, i)) / 0.01;
        end

    end

    for j = 1:size_Psi_2D_sp(2)

        for i = 1:size_Psi_2D_sp(1) - 1
            dpsi_dz(:, i) = (Psi_2D_sp(i + 1, j) - Psi_2D_sp(i, j)) / 0.01;
        end

    end

    size(dpsi_dr)
    size(dpsi_dz)
