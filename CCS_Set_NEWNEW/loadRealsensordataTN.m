%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS] = loadRealsensordata(PARAM)

    sensordata_B0 = fileread([PARAM.input_file_directory '/Sensor_B.txt']);
    sensordata_B = strsplit(sensordata_B0, {'\n', '\t', '\r'});
    sensornum_B = (length(sensordata_B) - 1) / 5 - 1;

    % 容器の外のセンサーを含めるなら1(=true)
    flag = 1;

    chnum = 0;
    countTPRB = 0;
    countNPRB = 0;

    for i = 1:sensornum_B
        R = str2double(sensordata_B{6 + (i - 1) * 5});
        Z = str2double(sensordata_B{7 + (i - 1) * 5});
        %PSI = str2double(sensordata_B{8+(i-1)*5});
        BZ = str2double(sensordata_B{9 + (i - 1) * 5});
        BR = str2double(sensordata_B{10 + (i - 1) * 5});
        % BR = 0;

        % 容器の内側か判定
        if (flag || R < 0.67)
            chnum = chnum + 1;

            % 接線方向・法線方向の判定
            if (abs(Z) < 0.95)
                countTPRB = countTPRB + 1;
                SENSOR_TPRB.R(countTPRB) = R;
                SENSOR_TPRB.Z(countTPRB) = Z;
                SENSOR_TPRB.TET(countTPRB) = atan2(BZ, BR);
                SENSOR_TPRB.TPRB(countTPRB) = sqrt(BR^2 + BZ^2);
                % SENSOR_TPRB.TPRB(countTPRB) = BZ;
                SENSOR_TPRB.ITYPE(countTPRB) = 1; % 使っていない
            else
                countNPRB = countNPRB + 1;
                SENSOR_NPRB.R(countNPRB) = R;
                SENSOR_NPRB.Z(countNPRB) = Z;
                SENSOR_NPRB.TET(countNPRB) = atan2(BZ, BR);
                SENSOR_NPRB.NPRB(countNPRB) = sqrt(BR^2 + BZ^2);
                % SENSOR_NPRB.NPRB(countNPRB) = BR;
                SENSOR_NPRB.ITYPE(countNPRB) = 1; % 使っていない
            end

            SENSOR_TPRB.TNPRB(chnum) = sqrt(BR^2 + BZ^2);
            RR(chnum) = R;
            ZZ(chnum) = Z;
            B(chnum) = sqrt(BR^2 + BZ^2);

        end

    end

    SENSOR_TPRB.NUM = countTPRB;
    SENSOR_NPRB.NUM = countNPRB;
    disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);
    disp(['Number of NPRB =  ' num2str(SENSOR_NPRB.NUM)]);

    % インボード側の磁束のプロットと、CCS面の決定2021年4月18日
    % スプライン補完でCCS面決定2021/05/06
    figure()
    RR0index = RR < 0.15;
    B = B(RR0index); ZZ = ZZ(RR0index);
    plot(ZZ, B);
    ZZ = ZZ'; B = B';
    f = fit(ZZ, B, 'smoothingspline', 'SmoothingParam', 1);
    hold on
    fnplt(f.p)
    x = fnzeros(fnder(f.p));
    x = unique(x(:));
    y = fnval(f.p, x);
    plot(x, y, 'o')
    xlim([-1 1]);
    matrix = [x y];
    matrix = sortrows(matrix, 2, 'descend');
    plot(matrix(1, 1), matrix(1, 2), "*");
    plot(matrix(2, 1), matrix(2, 2), "*");
    hold off
    lmax = islocalmax(B);
    view(90, 90);

    if length(ZZ(lmax)) > 1

        for i = 1:2
            CCS(i) = matrix(i, 1);
        end

    else
        size(matrix)
        CCS(1) = matrix(1, 1);
    end

    % ここまで

    sensordata_Flux0 = fileread([PARAM.input_file_directory '/Sensor_Flux.txt']);
    sensordata_Flux = strsplit(sensordata_Flux0, {'\n', '\t', '\r'});
    sensornum_Flux = (length(sensordata_Flux) - 1) / 5 - 1;

    chnum = 0;

    for i = 1:sensornum_Flux
        R = str2double(sensordata_Flux{6 + (i - 1) * 5});
        Z = str2double(sensordata_Flux{7 + (i - 1) * 5});
        PSI = str2double(sensordata_Flux{8 + (i - 1) * 5});
        BZ = str2double(sensordata_Flux{9 + (i - 1) * 5});
        BR = str2double(sensordata_Flux{10 + (i - 1) * 5});

        if (flag || R < 0.67)

            chnum = chnum + 1;

            SENSOR_FLXLP.R(chnum) = R;
            SENSOR_FLXLP.Z(chnum) = Z;
            SENSOR_FLXLP.FLXLP(chnum) = PSI;
            SENSOR_FLXLP.TET(chnum) = 0.0D0;
            SENSOR_FLXLP.ITYPE(chnum) = 0;
        end

    end

    SENSOR_FLXLP.NUM = chnum;
    disp(['Number of FLXLP = ' num2str(SENSOR_FLXLP.NUM)]);

end

%% loadsensordata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
