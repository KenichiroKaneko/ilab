function CCS_UTST_matlab(inputfile)
    %% *****************************************************************
    %% *  CAUCHY CONDITION SURFACE METHOD  (quadratic element version) *
    %% *****************************************************************

    format long e;
    close all;

    %MXCCS = 20;          % MAX! / NUMBER OF ELEMENTS ON THE CCS
    %MXINT = 10000;       % MAX! / NUMBER OF INTERNAL POINTS

    PARAM = loadinputfile(inputfile);

    REF = loadreference(PARAM);

    % 2021/05/06実際のセンサー配置にするかどうか
    % [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS] = loadsensordata(PARAM);
    [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS] = loadRealsensordata(PARAM);
    % 2021/05/06

    % CCS面の自動決定
    for i = 1:PARAM.CCS
        PARAM.Z0(i) = CCS(i);
    end

    WALL = loadwalldata(PARAM);
    % 2021/05/08
    % save('WALL',  'WALL');
    % error('Wall saved');
    % 2021/05/08

    ExtCOIL = loadcoildata(PARAM);

    FFDAT = makeFFdata(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP);

    GHR = zeros(1, 300);
    GHZ = zeros(1, 300);
    CCSDAT = makeCCSdata(PARAM, GHR, GHZ);

    FF = FFDAT;
    % 2021/05/07
    % save("FF", "FF");
    % error('FFsaved');
    % 2021/05/07
    FF(end + 1:end + sum(CCSDAT.NCCN)) = 0.0;
    if PARAM.IPCONST
        FF(end + 1) = -66410 * 4.0 * pi * 1.0e-7;
    end

    % 2021/05/08
    save('vars', 'WALL', 'FF', 'CCSDAT', 'PARAM', 'CCS', 'SENSOR_NPRB', 'SENSOR_TPRB', 'SENSOR_FLXLP');
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
    save('vars_afterFF', 'AA', 'FF', 'PSIFLX', 'FC', 'BZ', 'BR');
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
        SVD_MT_matlab2(PARAM, AA, FF, FC, 0, 0.0D0, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, FLXLP);

    half_norm = sqrt((sum((FFOUT).^2)));
    fprintf('%s%d\r\n', 'norm of the solution vector = ', half_norm);


    EDDYP(FFOUT, PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

    % plot eddy current
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
    JEDDY = dlmread(strcat([PARAM.input_file_directory '/jeddy.txt']));
    plot(JEDDY(:, 1), JEDDY(:, 2), '-b')

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
    [PSI, DELGE, RCCS, ZCCS, XPSI] = INTER(PARAM, 1, GETA, CR, CZ, FFOUT, ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, NINT);

    CCR = unique(CR);
    CCZ = unique(CZ);
    CCR(1) = [];
    psi = reshape(PSI(1:numel(CCR) * numel(CCZ)), numel(CCZ), numel(CCR));
    figure
    contour(CCR, CCZ, psi, '-k', 'LevelStep', 0.0003);
    hold on
    contour(REF.R, REF.Z, REF.Flux, '--m', 'LevelStep', 0.0003); % ????
    % plot(CCSDAT.RGI, CCSDAT.ZGI);
    hold off
    xlabel({'r (m)'});
    ylabel({'z (m)'});
    title("M-CCS")
    axis equal

    %% Data positions for 2D plot
    PLOT.R = 0.1:0.01:0.9;
    PLOT.Z = -1:0.02:1;

    fclose('all');
end

%% Main kokomade
