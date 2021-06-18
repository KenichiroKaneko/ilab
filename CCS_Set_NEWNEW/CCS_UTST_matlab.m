function CCS_UTST_matlab(inputfile)
    %% *****************************************************************
    %% *  CAUCHY CONDITION SURFACE METHOD  (quadratic element version) *
    %% *****************************************************************

    format long e;
    close all;

    %MXCCS = 20;          % MAX! / NUMBER OF ELEMENTS ON THE CCS
    %MXINT = 10000;       % MAX! / NUMBER OF INTERNAL POINTS

    PARAM = loadinputfile(inputfile);

    CONFIG = loadConfigure("CCS_input/config1.dat");

    [REF, JEDDY] = loadreference(PARAM, CONFIG);

    [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R] = loadSensorDataX(PARAM, CONFIG);

    for i = 1:PARAM.CCS
        PARAM.Z0(i) = CCS_Z(i);
        PARAM.R0(i) = CCS_R(i);
    end

    WALL = loadwalldata(PARAM);

    ExtCOIL = loadcoildata(PARAM);

    FFDAT = makeFFdata(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP);
    save("vars_aftermakeFFdata");

    GHR = zeros(1, 300);
    GHZ = zeros(1, 300);
    CCSDAT = makeCCSdata(PARAM, GHR, GHZ);

    % 各センサーポジションを表示する
    if CONFIG.ShowFig
        dispSensorPosition(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, REF);
    end

    FF = FFDAT;
    FF(end + 1:end + sum(CCSDAT.NCCN)) = 0.0;

    if PARAM.IPCONST
        FF(end + 1) = -66410 * 4.0 * pi * 1.0e-7;
    end

    % 2021/05/08
    save("vars_afterLoadData");
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
        FORM(PARAM, CONFIG, AA, FF, ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

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

    [C, W, U, V, FFOUT, XBFR, XMT] = ...
        SVD_MT_matlab2(PARAM, CONFIG, AA, FF, FC, 0, 0.0D0, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, FLXLP);

    % save("vars_afterSVD");
    % error('error description')

    half_norm = sqrt((sum((FFOUT).^2)));
    fprintf('%s%d\r\n', 'norm of the solution vector = ', half_norm);

    % EDDYP(FFOUT, PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % plot eddy current
    if CONFIG.ShowFig
        DISF = dlmread([PARAM.output_file_directory '/EddyCurrentProfile.txt']);
        figure('Name', 'Eddy Current Plofile', 'NumberTitle', 'off')
        plot(DISF(:, 1), DISF(:, 2), '-ko', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 2)
        xlabel({'Distance (m)'});
        ylabel({'Eddy Current Density (MA/m)'});
        axis([0 DISF(end, 1) -0.2 0.2])
        set(gca, 'FontSize', 14);

        axis([0 DISF(end, 1) -0.4 0.4])
        set(gca, 'FontSize', 14);
        hold on
        plot(JEDDY(:, 1), JEDDY(:, 2), '-b')
    end

    if 0
        % plot sensor position
        VV = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);
        SEN0 = dlmread([PARAM.temporary_file_directory '/SENPOS0.txt']);
        SEN1 = dlmread([PARAM.temporary_file_directory '/SENPOS1.txt']);

        figure('Name', 'Sensor Position', 'NumberTitle', 'off')
        subplot(1, 2, 1);
        plot(VV(:, 1), VV(:, 2), '-k');
        hold on
        plot(SEN0(:, 1), SEN0(:, 2), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4)
        plot(SEN1(:, 1), SEN1(:, 2), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4)
        title('Sensor Position')
        xlabel('r [m]');
        ylabel('z [m]');
        axis equal
        axis([0 1 -1.2 1.2])
        set(gca, 'FontSize', 14);

        % plot CV segment
        %COIL = dlmread('output\@UTST_CoilGeom.txt');
        VVMESH = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);
        VVNODE = dlmread([PARAM.temporary_file_directory '/VacuumVesselNodePoints.txt']);
        VVSEG = dlmread([PARAM.temporary_file_directory '/VacuumVesselSegments.txt']);
        subplot(1, 2, 2);
        %figure('Name','CVsegment','NumberTitle','off')
        plot(VV(:, 1), VV(:, 2), '-k')
        hold on
        plot(VVMESH(:, 1), VVMESH(:, 2), 'bx', 'MarkerSize', 8)
        hold on
        plot(VVNODE(:, 1), VVNODE(:, 2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
        % plot(COIL(:,1),COIL(:,2),'s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
        % plot(RES(1:3),ZES(1:3),'ko',RES(1:3),ZES(1:3),'-k')
        title('CVsegment');
        xlabel('r [m]');
        ylabel('z [m]');
        axis equal
        axis([0 1 -1.2 1.2])
        set(gca, 'FontSize', 14);

        for i = 1:PARAM.CCS
            RCCSSP(i, :) = spline(1:CCSDAT.NCCN + 1, horzcat(CCSDAT.RCCN(i, 1:CCSDAT.NCCN), CCSDAT.RCCN(i, 1)), 1:1/10:CCSDAT.NCCN + 1);
            ZCCSSP(i, :) = spline(1:CCSDAT.NCCN + 1, horzcat(CCSDAT.ZCCN(i, 1:CCSDAT.NCCN), CCSDAT.ZCCN(i, 1)), 1:1/10:CCSDAT.NCCN + 1);

            plot(RCCSSP(i, :), ZCCSSP(i, :), '-k')
            plot(CCSDAT.RCCS(i, 1:CCSDAT.NCCS), CCSDAT.ZCCS(i, 1:CCSDAT.NCCS), 'bx', 'MarkerSize', 8)
            plot(CCSDAT.RCCN(i, 1:CCSDAT.NCCN), CCSDAT.ZCCN(i, 1:CCSDAT.NCCN), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
        end

    end

    MINR = 10;
    MAXR = 90;
    MINZ = -100;
    MAXZ = 100;
    ICRE = 1;
    JCRE = 2;
    NINT = 0;

    % 再構成する磁束のサイズ
    % cm単位の範囲の指定
    % MINR = 10.815;
    % MAXR = 88.80;
    % MINZ = -99.85;
    % MAXZ = 99.85;
    % % メッシュの間隔delr,delz
    % ICRE = 9.747920133111480e-02;
    % JCRE = 9.827755905511811e-02;
    % NINT = 0;
    % Nz = 2033;
    % Nr = 800;
    % zmin = -9.985000000000001e-01;
    % zmax = 9.985000000000001e-01;
    % rmin = 1.081500000000000e-01;
    % rmax = 0.8880;
    % delr = 9.747920133111480e-04;
    % delz = 9.827755905511811e-04;

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

    if CONFIG.ShowFig
        figure
        contour(CCR, CCZ, psi, '-k', 'LevelStep', 0.0003);
        hold on
        contour(REF.R, REF.Z, REF.Flux, '--m', 'LevelStep', 0.0003); % ????
        % scatter(CCSDAT.RCCN, CCSDAT.ZCCN, 'o')
        plot(CCSDAT.RGI, CCSDAT.ZGI);
        hold off
        xlabel({'r (m)'});
        ylabel({'z (m)'});
        title("Reconstructed flux")
        axis equal

        %% Data positions for 2D plot
        PLOT.R = 0.1:0.01:0.9;
        PLOT.Z = -1:0.02:1;
    end

    % 評価関数
    MSE = EVALUATE(psi, REF, PARAM, CONFIG, CCSDAT, CCR, CCZ)
    % save("vars_result_" + PARAM.input_file_directory, "psi", "REF", "PARAM", "CCR", "CCZ", "CCSDAT", "MSE");

    fclose('all');
end

%% Main kokomade
