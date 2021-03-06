%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R] = loadsensordata(PARAM)

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

    % 2021/5/24 SVD_MT_matlab2での整合性
    SENSOR_TPRB.TNPRB = SENSOR_TPRB.TPRB;

    SENSOR_TPRB.NUM = chnum * 2;
    disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);

    for i = 1:PARAM.CCS
        CCS_Z(i) = PARAM.Z0(i);
        CCS_R(i) = PARAM.R0(i);
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

    % save('vars_RZ')

    % error('error description', A1)
    %
    %     %%  *************************************************************************
    %     %%     Generation of CCS input data
    %     %%  *************************************************************************
end

%% loadsensordata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% loadRealsensordataをここに保存しておく。2021/10/05
% function [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R] = loadRealsensordata(PARAM, CONFIG)

%     sensordata_B0 = fileread([PARAM.input_file_directory '/Sensor_B.txt']);
%     sensordata_B = strsplit(sensordata_B0, {'\n', '\t', '\r'});
%     sensornum_B = (length(sensordata_B) - 1) / 5 - 1;

%     % 容器の外のセンサーを含めるなら1(=true)
%     flag = 1;
%     chnum = 0;

%     % 磁場センサーを読み込む
%     for i = 1:sensornum_B
%         R = str2double(sensordata_B{6 + (i - 1) * 5});
%         Z = str2double(sensordata_B{7 + (i - 1) * 5});
%         %PSI = str2double(sensordata_B{8+(i-1)*5});
%         BZ = str2double(sensordata_B{9 + (i - 1) * 5});
%         BR = str2double(sensordata_B{10 + (i - 1) * 5});

%         % 容器の内側か判定している
%         if (flag || R < 0.67)
%             chnum = chnum + 1;
%             SENSOR_TPRB.R(chnum) = R;
%             SENSOR_TPRB.Z(chnum) = Z;
%             SENSOR_TPRB.TET(chnum) = atan2(BZ, BR);
%             SENSOR_TPRB.TPRB(chnum) = sqrt(BR^2 + BZ^2);
%             SENSOR_TPRB.ITYPE(chnum) = 1; % 使っていない
%             RR(chnum) = R;
%             ZZ(chnum) = Z;
%             B(chnum) = sqrt(BR^2 + BZ^2);
%             SENSOR_TPRB.TNPRB(chnum) = sqrt(BR^2 + BZ^2);
%         end

%     end

%     SENSOR_TPRB.NUM = chnum;
%     disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);

%     % インボード側の磁束のプロットと、CCS面の決定2021年4月18日
%     % スプライン補完でCCS面決定2021/05/06
%     if CONFIG.DetermineCCSZPos
%         figure()
%         RR0index = RR < 0.15;
%         B = B(RR0index); ZZ = ZZ(RR0index);
%         plot(ZZ, B);
%         ZZ = ZZ'; B = B';
%         f = fit(ZZ, B, 'smoothingspline', 'SmoothingParam', 1);
%         hold on
%         fnplt(f.p)
%         x = fnzeros(fnder(f.p)); x = unique(x(:));
%         y = fnval(f.p, x);
%         plot(x, y, 'o')
%         xlim([-1 1]); ylabel("ψ[Wb]"); xlabel("Z[m]");
%         matrix = [x y];
%         matrix = sortrows(matrix, 2, 'descend');
%         plot(matrix(1, 1), matrix(1, 2), "*"); plot(matrix(2, 1), matrix(2, 2), "*");
%         hold off
%         lmax = islocalmax(B);
%         view(90, 90);

%         if length(ZZ(lmax)) > 1

%             for i = 1:2
%                 CCS_Z(i) = matrix(i, 1);
%             end

%         else
%             size(matrix)
%             CCS_Z(1) = matrix(1, 1);
%         end

%     else

%         for i = 1:PARAM.CCS
%             CCS_Z(i) = PARAM.Z0(i);
%             CCS_R(i) = PARAM.R0(i);
%         end

%     end

%     % CCS面の自動決定R方向
%     if CONFIG.DetermineCCSRPos
%         CCS_R = CalcPlasmaCenter(PARAM, CONFIG, CCS_Z);
%         % CCS_R = CCS_R * 0.9;
%         % CCS_R = 0.295
%     end

%     % ここまで

%     %% No NPRB
%     SENSOR_NPRB.NUM = 0;
%     SENSOR_NPRB.R = [];
%     SENSOR_NPRB.Z = [];
%     SENSOR_NPRB.TET = [];
%     SENSOR_NPRB.NPRB = [];
%     SENSOR_NPRB.ITYPE = [];

%     sensordata_Flux0 = fileread([PARAM.input_file_directory '/Sensor_Flux.txt']);
%     sensordata_Flux = strsplit(sensordata_Flux0, {'\n', '\t', '\r'});
%     sensornum_Flux = (length(sensordata_Flux) - 1) / 5 - 1;

%     chnum = 0;

%     % 磁束データを読み込む
%     for i = 1:sensornum_Flux
%         R = str2double(sensordata_Flux{6 + (i - 1) * 5});
%         Z = str2double(sensordata_Flux{7 + (i - 1) * 5});
%         PSI = str2double(sensordata_Flux{8 + (i - 1) * 5});
%         BZ = str2double(sensordata_Flux{9 + (i - 1) * 5});
%         BR = str2double(sensordata_Flux{10 + (i - 1) * 5});

%         if (flag || R < 0.67)

%             chnum = chnum + 1;

%             SENSOR_FLXLP.R(chnum) = R;
%             SENSOR_FLXLP.Z(chnum) = Z;
%             % SENSOR_FLXLP.FLXLP(chnum) = PSI;
%             SENSOR_FLXLP.FLXLP(chnum) = PSI / (pi * R^2); % 21/05/30
%             SENSOR_FLXLP.TET(chnum) = 0.0D0;
%             SENSOR_FLXLP.ITYPE(chnum) = 0;
%         end

%     end

%     SENSOR_FLXLP.NUM = chnum;
%     disp(['Number of FLXLP = ' num2str(SENSOR_FLXLP.NUM)]);

% end
