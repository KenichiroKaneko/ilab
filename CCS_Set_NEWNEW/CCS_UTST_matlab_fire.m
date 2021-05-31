function CCS_UTST_matlab(inputfile)
    %% *****************************************************************
    %% *  CAUCHY CONDITION SURFACE METHOD  (quadratic element version) *
    %% *****************************************************************

    format long e;
    close all;

    %MXCCS = 20;          % MAX! / NUMBER OF ELEMENTS ON THE CCS
    %MXINT = 10000;       % MAX! / NUMBER OF INTERNAL POINTS

    PARAM = loadinputfile(inputfile);
    PARAM.dispFigures = 0;

    REF = loadreference(PARAM);

    % 2021/05/06実際のセンサー配置にするかどうか
    % オリジナルのコード、入力は、末尾にRがないもの、UTST_numel_5,9が対応
    % [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z] = loadsensordata(PARAM);

    % 実際のセンサー位置で再構成するコード、入力は末尾がRと、実験データが対応(例)UTST_numel_2033R、180515_010_t9650
    [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R] = loadRealsensordata(PARAM);

    % 法線方向、接線方向を考慮して読み込むためのコード、入力は末尾がRと、実験データが対応
    % [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z,CCS_R] = loadRealsensordataTN(PARAM);

    for i = 1:PARAM.CCS
        PARAM.Z0(i) = CCS_Z(i);
        PARAM.R0(i) = CCS_R(i);
    end

    WALL = loadwalldata(PARAM);

    ExtCOIL = loadcoildata(PARAM);

    FFDAT = makeFFdata(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP);
    % save("vars_aftermakeFFdata");

    GHR = zeros(1, 300);
    GHZ = zeros(1, 300);
    CCSDAT = makeCCSdata(PARAM, GHR, GHZ);

    % 各センサーポジションを表示する
    if PARAM.dispFigures
        dispSensorPosition(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, REF);
    end

    FF = FFDAT;
    FF(end + 1:end + sum(CCSDAT.NCCN)) = 0.0;

    if PARAM.IPCONST
        FF(end + 1) = -66410 * 4.0 * pi * 1.0e-7;
    end

    % 2021/05/08
    % save("vars_afterLoadData");
    % error('vars saved');
    % 2021/05/08
    %% *************************************************************
    %%      iteration Start
    DELGE = 0.0;
    OLDEL = 1.0;
    GETA = 0.0; %! 下駄の初期値をゼロにする
    OLDGT = GETA;
    OMGA = 1.9;

    if (PARAM.ITSKP > 0)
        ITMX = 1;
        fprintf('下駄サーチを (1) SVD_MT内部のみでやります\n');
    else
        fprintf('下駄サーチを (0) INTERの反復でもやります\n');
    end

    if (PARAM.IPCONST == 1)
        NMAX = SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN) + 1; % ushikiip
    else
        NMAX = SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN); % ushikiip
    end

    JMAX = sum(CCSDAT.NCCN) + sum(CCSDAT.NCCN) + sum(WALL.KNN) + sum(WALL.KSN);

    AA = zeros(NMAX, JMAX);
    [FC, BR, BZ, PSIFLX, PSIC, AA, FF] = ...
        FORM(PARAM, AA, FF, ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % 2021/05/08
    % save('vars_afterForm');
    % error('saved vars')
    % 2021/05/08

    %% INOMOTO start
    fluxfactor = 10;

    FF(SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM) = FF(SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM) * fluxfactor;
    AA(SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM, :) = AA(SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM, :) * fluxfactor;
    FC(SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM) = FC(SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM) * fluxfactor;

    FLXLP = FFDAT(SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM) * fluxfactor;

    %% INOMOTO end

    %% Solve Matrix Equation
    %% Singular Value Decompoition

    fprintf('NCCN/KNN/KSN = %d %d %d\n', sum(CCSDAT.NCCN), sum(WALL.KNN), sum(WALL.KSN));

    fireNote = fopen([PARAM.temporary_file_directory '/@fireNote_' PARAM.input_file_directory '.txt'], 'w');
    % for i = 1:length(FF)
    for i = 1:NMAX
        superKUP = i;

        try
            [C, W, U, V, FFOUT, XBFR, XMT] = ...
                SVD_MT_matlab_fire(PARAM, AA, FF, FC, 0, 0.0D0, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, FLXLP, superKUP);

            % % save("vars_afterSVD");
            % error('error description')

            EDDYP(FFOUT, PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

            MINR = 10;
            MAXR = 90;
            MINZ = -100;
            MAXZ = 100;
            ICRE = 1;
            JCRE = 2;
            NINT = 0;

            for I = MINR:ICRE:MAXR
                NCOUNT = 0;
                CCR = I / 100.0;

                for J = MINZ:JCRE:MAXZ
                    NINT = NINT + 1;
                    NCOUNT = NCOUNT + 1;
                    CR(NINT) = CCR;
                    CZ(NINT) = J / 100.0;
                end

            end

            % 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
            [PSI, DELGE, RCCS, ZCCS, XPSI] = INTER(PARAM, 1, GETA, CR, CZ, FFOUT, ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, NINT);

            CCR = unique(CR);
            CCZ = unique(CZ);
            CCR(1) = [];
            psi = reshape(PSI(1:numel(CCR) * numel(CCZ)), numel(CCZ), numel(CCR));

            % 評価関数
            MSE(i) = EVALUATE(psi, REF, PARAM, CCSDAT, CCR, CCZ);
            % save("vars_result_" + PARAM.input_file_directory, "psi", "REF", "PARAM", "CCR", "CCZ", "CCSDAT", "MSE");
            fprintf(fireNote, '%3d   KUP=%3d   MSE=%.3e\n', i, superKUP, MSE(i));
            fprintf('%3d   KUP=%3d   MSE=%.3e\n', i, superKUP, MSE(i));
        catch ME
            fprintf('error at i = %d', i);
            disp(ME)
            MSE(i) = -1;
            fprintf(fireNote, '%3d   KUP=%3d   MSE=%.3e\n', i, superKUP, MSE(i));
            fprintf('%3d   KUP=%3d   MSE=%.3e\n', i, superKUP, MSE(i));
        end

    end

    save("vars_result_" + PARAM.input_file_directory, "psi", "REF", "PARAM", "CCR", "CCZ", "CCSDAT", "MSE")

    fclose('all');
end

%% Main kokomade
