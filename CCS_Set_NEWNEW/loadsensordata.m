%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS] = loadsensordata(PARAM)

    sensordata_B0 = fileread([PARAM.input_file_directory '/Sensor_B.txt']);
    sensordata_B = strsplit(sensordata_B0, {'\n', '\t', '\r'});
    sensornum_B = (length(sensordata_B) - 1) / 5 - 1;

    chnum = 0;

    for i = 1:sensornum_B
        R = str2double(sensordata_B{6 + (i - 1) * 5});
        Z = str2double(sensordata_B{7 + (i - 1) * 5});
        %PSI = str2double(sensordata_B{8+(i-1)*5});
        BZ = str2double(sensordata_B{9 + (i - 1) * 5});
        BR = str2double(sensordata_B{10 + (i - 1) * 5});

        if ((R < 0.57 + 0.0015 && R > 0.27 - 0.0015) || (Z < 0.239 && R > 0.67) || (R < 0.689 && R > 0.629 - 0.0015))
        else
            chnum = chnum + 1;

            SENSOR_TPRB.R(chnum * 2 - 1) = R;
            SENSOR_TPRB.R(chnum * 2) = R;
            SENSOR_TPRB.Z(chnum * 2 - 1) = Z;
            SENSOR_TPRB.Z(chnum * 2) = -Z;
            SENSOR_TPRB.TET(chnum * 2 - 1) = atan2(BZ, BR);
            SENSOR_TPRB.TET(chnum * 2) = atan2(BZ, -BR);
            SENSOR_TPRB.TPRB(chnum * 2 - 1) = sqrt(BR^2 + BZ^2);
            SENSOR_TPRB.TPRB(chnum * 2) = sqrt(BR^2 + BZ^2);
            SENSOR_TPRB.ITYPE(chnum * 2 - 1) = 1;
            SENSOR_TPRB.ITYPE(chnum * 2) = 1;
            RR(chnum) = R;
            ZZ(chnum) = Z;
            B(chnum) = sqrt(BR^2 + BZ^2);
        end

    end

    SENSOR_TPRB.NUM = chnum * 2;
    disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);

    % 2021/05/21
    figure()
    % RR = [RR RR];
    % ZZ = [ZZ -ZZ];
    scatter(RR, ZZ)
    figure()
    hold on
    % B = [B B];
    RR0index = RR < 0.15;
    B = B(RR0index); ZZ = ZZ(RR0index);
    B = [B(1) * 2 - B(2) B B(end - 1) * 2 - B(end)]
    ZZ = [ZZ(1) * 2 - ZZ(2) ZZ ZZ(end) * 2 - ZZ(end - 1)]
    ZZ = ZZ'; B = B';
    plot(fliplr(ZZ), B);
    f = fit(ZZ, B, 'smoothingspline', 'SmoothingParam', 1);
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
    lmax = islocalmax(B)
    lmax(1) = [];
    lmax(end) = [];
    lmax

    view(90, 90);

    CCS1 = matrix(1, 1);
    CCS2 = matrix(2, 1);

    if length(ZZ(lmax)) > 1

        CCS(1) = max(CCS1, CCS2);
        CCS(2) = min(CCS1, CCS2);

    else
        size(matrix)
        CCS(1) = matrix(1, 1);
    end

    % 2021/05/21

    %% No NPRB
    SENSOR_NPRB.NUM = 0;
    SENSOR_NPRB.R = [];
    SENSOR_NPRB.Z = [];
    SENSOR_NPRB.TET = [];
    SENSOR_NPRB.NPRB = [];
    SENSOR_NPRB.ITYPE = [];

    sensordata_Flux0 = fileread([PARAM.input_file_directory '/Sensor_Flux.txt']);
    sensordata_Flux = strsplit(sensordata_Flux0, {'\n', '\t', '\r'});
    sensornum_Flux = (length(sensordata_Flux) - 1) / 5 - 1;

    chnum = 0;

    for i = 1:sensornum_Flux
        R = str2double(sensordata_Flux{6 + (i - 1) * 5});
        Z = str2double(sensordata_Flux{7 + (i - 1) * 5});
        PSI = str2double(sensordata_Flux{8 + (i - 1) * 5});
        %BZ = str2double(sensordata_Flux{8+(i-1)*5});
        %BR = str2double(sensordata_Flux{10+(i-1)*5});

        chnum = chnum + 1;

        SENSOR_FLXLP.R(chnum * 2 - 1) = R;
        SENSOR_FLXLP.R(chnum * 2) = R;
        SENSOR_FLXLP.Z(chnum * 2 - 1) = Z;
        SENSOR_FLXLP.Z(chnum * 2) = -Z;
        SENSOR_FLXLP.FLXLP(chnum * 2 - 1) = PSI;
        SENSOR_FLXLP.FLXLP(chnum * 2) = PSI;
        SENSOR_FLXLP.TET(chnum * 2 - 1) = 0.0D0;
        SENSOR_FLXLP.TET(chnum * 2) = 0.0D0;
        SENSOR_FLXLP.ITYPE(chnum * 2 - 1) = 0;
        SENSOR_FLXLP.ITYPE(chnum * 2) = 0;
    end

    SENSOR_FLXLP.NUM = chnum * 2;
    disp(['Number of FLXLP = ' num2str(SENSOR_FLXLP.NUM)]);

    if PARAM.IUTST == 5
        SENSOR_FLXLP.FLXLP = SENSOR_FLXLP.FLXLP - 0.0042459;
    end

    save('vars_RZ')
    % error('error description', A1)
    %
    %     %%  *************************************************************************
    %     %%     Generation of CCS input data
    %     %%  *************************************************************************
end

%% loadsensordata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
