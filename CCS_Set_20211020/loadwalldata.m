%% loadwalldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadwalldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadwalldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadwalldata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function WALL = loadwalldata(PARAM)

    %% Here, the geometory should be derived from input file
    %     R = [0.694D0, 0.694D0, 0.5985D0, 0.5985D0, 0.10815D0, 0.10815D0,...
    %         0.10815D0, 0.5985D0, 0.5985D0, 0.694D0,  0.694D0];
    %     Z = [0.0D0, 0.285D0, 0.285D0, 0.9985D0, 0.9985D0, 0.0D0,...
    %         -0.9985D0, -0.9985D0, -0.285D0, -0.285D0 , 0.0D0];

    %% CONDSHL_UTST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % 導体上の渦電流節点位置データの設定と読み取り
    %function [KNE,KNN,REV,ZEV,RES,ZES,KSE,KSN,REVN,ZEVN] = CONDSHL_UTST(WAHAHA,Nedp,NONC)

    MAXM = 13;
    Z = 0.15;
    INBOAD = 0.5;

    MSEC = [1, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 1];

    RSEC(1) = 0.694; ZSEC(1) = 0.0;
    RSEC(2) = 0.694; ZSEC(2) = 0.285;
    RSEC(3) = 0.5985; ZSEC(3) = 0.285;
    RSEC(4) = 0.5985; ZSEC(4) = 0.285 + Z;
    RSEC(5) = 0.5985; ZSEC(5) = 0.9985;
    RSEC(6) = 0.10815; ZSEC(6) = 0.9985;
    RSEC(7) = 0.10815; ZSEC(7) = INBOAD;
    RSEC(8) = 0.10815; ZSEC(8) = -INBOAD;
    RSEC(9) = 0.10815; ZSEC(9) = -0.9985;
    RSEC(10) = 0.5985; ZSEC(10) = -0.9985;
    RSEC(11) = 0.5985; ZSEC(11) = -0.285 - Z;
    RSEC(12) = 0.5985; ZSEC(12) = -0.285;
    RSEC(13) = 0.694; ZSEC(13) = -0.285;
    RSEC(14) = 0.694; ZSEC(14) = 0.0;
    %% Modified_MI 20211103
    RSEC(14) = []; ZSEC(14) = [];
    RSEC(1) = []; ZSEC(1) = [];
    % MSEC = [1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1]; %original
    MSEC = [1, 1, 2, 2, 2, 3, 2, 2, 2, 1, 1]; % ほどほどにノード点を増やした
    MAXM = 11;
    %% Modified_MI 20211103 kokomade
    % % kaneko
    % RSEC_list = [0.7, 0.6, 0.6, 0.1];
    % RSEC_list = [RSEC_list fliplr(RSEC_list)];
    % ZSEC_list = [0.285, 0.285, 0.9985, 0.9985];
    % ZSEC_list = [ZSEC_list -fliplr(ZSEC_list)];

    % MAXM = length(RSEC_list) - 1;

    % for i = 1:MAXM + 1
    %     RSEC(i) = RSEC_list(i);
    %     ZSEC(i) = ZSEC_list(i);
    % end

    % MSEC = ones(MAXM);
    % % MSEC = [1, 5, 1, 1, 1, 5, 1];
    % MSEC = MSEC + 1;
    % kaneko
    fp = fopen([PARAM.temporary_file_directory '/VacuumVesselSegments.txt'], 'w');

    for i = 1:MAXM + 1
        fprintf(fp, '%d %d\n', RSEC(i), ZSEC(i));
    end

    fclose(fp);

    % 真空容器上の渦電流節点の設定
    % KNE=真空容器の分割境界要素数,  KNM=真空容器上のメッシュ点数,  KNN=真空容器上の節点数

    II = 1;
    WALL.REV(II) = RSEC(1);
    WALL.ZEV(II) = ZSEC(1);

    WALL.KNE = 0;

    for I = 1:MAXM
        WALL.KNE = WALL.KNE + MSEC(I);
        CM = 2 * MSEC(I);
        DELR = (RSEC(I + 1) - RSEC(I)) / CM;
        DELZ = (ZSEC(I + 1) - ZSEC(I)) / CM;

        for J = 1:2 * MSEC(I)
            II = II + 1;
            WALL.REV(II) = WALL.REV(II - 1) + DELR;
            WALL.ZEV(II) = WALL.ZEV(II - 1) + DELZ;
        end

    end

    KNM = II - 1;

    fp = fopen([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt'], 'w');

    for II = 1:KNM + 1
        fprintf(fp, '%d %d\n', WALL.REV(II), WALL.ZEV(II));
    end

    fclose(fp);

    %% ！！！！！！！　 非適合要素節点計算の予定地　　　！！！！！！！
    % 真空容器上の非適合要素節点座標の作成
    if (PARAM.NONC > 0)
        WALL.KNN = WALL.KNE * 3;
        I = 1:WALL.KNE;
        WALL.REVN(3 * I - 2) = (5.D0 * WALL.REV(2 * I - 1) + 5.D0 .* WALL.REV(2 * I) - WALL.REV(2 * I + 1)) / 9.D0;
        WALL.ZEVN(3 * I - 2) = (5.D0 * WALL.ZEV(2 * I - 1) + 5.D0 .* WALL.ZEV(2 * I) - WALL.ZEV(2 * I + 1)) / 9.D0;
        WALL.REVN(3 * I - 1) = WALL.REV(2 * I);
        WALL.ZEVN(3 * I - 1) = WALL.ZEV(2 * I);
        WALL.REVN(3 * I) = (5.D0 * WALL.REV(2 * I + 1) + 5.D0 .* WALL.REV(2 * I) - WALL.REV(2 * I - 1)) / 9.D0;
        WALL.ZEVN(3 * I) = (5.D0 * WALL.ZEV(2 * I + 1) + 5.D0 .* WALL.ZEV(2 * I) - WALL.ZEV(2 * I - 1)) / 9.D0;
    else
        WALL.KNN = KNM;
    end

    fp = fopen([PARAM.temporary_file_directory '/VacuumVesselNodePoints.txt'], 'w'); %113

    if (PARAM.NONC == 0)

        for II = 1:KNM + 1
            fprintf(fp, '%d %d\n', WALL.REV(II), WALL.ZEV(II));
        end

        fprintf('KNE=%d    KNM=%d\n', WALL.KNE, KNM);
        fprintf('%s\n', 'The end meshpoint agrees with the start meshpoint');
    else

        for II = 1:WALL.KNN
            fprintf(fp, '%d %d\n', WALL.REVN(II), WALL.ZEVN(II));
        end

        fprintf('KNE/KNM/KNN = %d  %d  %d\n', WALL.KNE, KNM, WALL.KNN);
    end

    fclose(fp);

    % 安定化板上の渦電流節点の設定
    % KSE=安定化板の分割境界要素数,  KSN=安定化板上の節点数
    %***
    %cd      KSE=1
    %prompt = '電流シート上に境界要素を配置する？ (Yes/No)=(1/0)\n';
    %KSE = input(prompt);

    WALL.KSE = 0;
    %***
    WALL.KSN = WALL.KSE * 2 + 1;

    if (WALL.KSE == 0)
        WALL.KSN = 0;
    end

    if (WALL.KSE > 0)
        fp = fopen([PARAM.temporary_file_directory '/StabilizerPoints.txt'], 'w'); % 114
        WALL.RES(1) = 0.25;
        WALL.RES(2) = 0.28;
        WALL.RES(3) = 0.31;
        WALL.ZES(1) = 0;
        WALL.ZES(2) = 0;
        WALL.ZES(3) = 0;

        for I = 1:WALL.KSN
            fprintf(fp, '%d %d\n', WALL.RES(I), WALL.ZES(I));
        end

        fclose(fp);
    else
        WALL.RES = [];
        WALL.ZES = [];
    end

end

%% loadwalldata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
