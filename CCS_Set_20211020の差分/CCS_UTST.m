% function CCS_UTST()
%% *****************************************************************
%% *  CAUCHY CONDITION SURFACE METHOD  (quadratic element version) *
%% *****************************************************************

format long e;
close all;

% CONFIG : 2πR^2で割るかや、プロットの表示、実験データか数値解かを設定するファイル
% CONFIG = loadConfigure("CCS_input/config1.dat");
CONFIG = loadConfigure("CCS_input/configForExp.dat");

% 数値解はDataType=0
if CONFIG.DataType == 0
    PARAM = loadinputfile('input_m2033.txt');
    [REF, JEDDY] = loadreference(PARAM, CONFIG);
    WALL = loadwalldata(PARAM);
    ExtCOIL = loadcoildata(PARAM);
    [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R, REF] = loadSampleData(PARAM, CONFIG);
elseif CONFIG.DataType == 1
    PARAM = loadinputfile('input_5.txt');
    REF = 0;
    WALL = loadwalldata(PARAM);
    % 211018 #004-#006 がデータが取れたショット
    PARAM.shotnum = '006';
    PARAM.foldername = 'utst/211018';
    PARAM.time_CCS = 9.553;
    [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, ExtCOIL, CCS_Z, CCS_R] = loadLizzieData(PARAM, CONFIG);
    % ExtCOIL.I = [0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
end

% 実験値のFLXLPとReferenceの磁場データ
% AAAA = SENSOR_FLXLP.FLXLP;
% BBBB = SENSOR_TPRB.TPRB;
% SENSOR_FLXLP.FLXLP = AAAA;
% SENSOR_TPRB.TPRB = BBBB;

% 死んでそうな磁場センサーを消す
i = [19:24, 30, 35:40];
SENSOR_TPRB.R(i) = [];
SENSOR_TPRB.Z(i) = [];
SENSOR_TPRB.TET(i) = [];
SENSOR_TPRB.TPRB(i) = [];
SENSOR_TPRB.TNPRB(i) = [];
SENSOR_TPRB.NUM = length(SENSOR_TPRB.TPRB);
% error('error description', A1)

%
for i = 1:PARAM.CCS
    PARAM.Z0(i) = CCS_Z(i);
    PARAM.R0(i) = CCS_R(i);
end

FFDAT = makeFFdata(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP);
GHR = zeros(1, 300);
GHZ = zeros(1, 300);
CCSDAT = makeCCSdata(PARAM, GHR, GHZ);

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
    FORM(PARAM, CONFIG, AA, FF, ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

%% INOMOTO start
% 規格化に使うそれぞれのセンサーのインデックス
range_b = 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM;
range_FLXLP = SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM;
range_CCS = SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM + 1:SENSOR_NPRB.NUM + SENSOR_TPRB.NUM + SENSOR_FLXLP.NUM + CCSDAT.NCCN;

% 平均が１になるように規格化→センサー信号が全て０付近の場合悪くなりそう
% bfactor = 1 / mean(FF(range_b));
% fluxfactor = 1 / mean(FF(range_FLXLP));
% ccsfactor = 1 / mean(FF(range_CCS));

% 絶対値の最大値が１になるように規格化
% bfactor = 1 / max(abs(FF(range_b)));
% fluxfactor = 1 / max(abs(FF(range_FLXLP)));
% ccsfactor = 1 / max(abs(FF(range_CCS)));

% ノルムで割る規格化
% bfactor = 1 / norm(FF(range_b), 1);
% fluxfactor = 1 / norm(FF(range_FLXLP), 1);
% ccsfactor = 1 / norm(FF(range_CCS), 1);

% 手動の規格化
% bfactor = 1;
% fluxfactor = 10;
% ccsfactor = 1;

%
FF(range_b) = FF(range_b) * bfactor;
AA(range_b, :) = AA(range_b, :) * bfactor;
FC(range_b) = FC(range_b) * bfactor;
FF(range_FLXLP) = FF(range_FLXLP) * fluxfactor;
AA(range_FLXLP, :) = AA(range_FLXLP, :) * fluxfactor;
FC(range_FLXLP) = FC(range_FLXLP) * fluxfactor;
FLXLP = FFDAT(range_FLXLP) * fluxfactor;
FF(range_CCS) = FF(range_CCS) * ccsfactor;
AA(range_CCS, :) = AA(range_CCS, :) * ccsfactor;
FC(range_CCS) = FC(range_CCS) * ccsfactor;

%% INOMOTO end

%% Solve Matrix Equation
%% Singular Value Decompoition

fprintf('NCCN/KNN/KSN = %d %d %d\n', sum(CCSDAT.NCCN), sum(WALL.KNN), sum(WALL.KSN));
% SVD_MT_matlab2でL-curve法を実装してある
% [C, W, U, V, FFOUT, XBFR, XMT] = ...
%     SVD_MT_matlab(PARAM, AA, FF, FC, 0, 0.0D0, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, FLXLP);
[C, W, U, V, FFOUT, XBFR, XMT] = ...
SVD_MT_matlab2(PARAM, CONFIG, AA, FF, FC, 0, 0.0D0, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, FLXLP, FF);

half_norm = sqrt((sum((FFOUT).^2)));
fprintf('%s%d\r\n', 'norm of the solution vector = ', half_norm);

EDDYP(FFOUT, PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);

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

    if CONFIG.DataType == 0
        hold on
        plot(JEDDY(:, 1), JEDDY(:, 2), '-b')
    end

end

if CONFIG.ShowFig
    % plot sensor position
    % VV : 容器壁
    VV = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);

    figure('Name', 'Sensor Position', 'NumberTitle', 'off')
    subplot(1, 3, 1);
    plot(VV(:, 1), VV(:, 2), '-k');
    hold on
    plot(SENSOR_FLXLP.R, SENSOR_FLXLP.Z, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 4)
    title('FluxLoop Position')
    xlabel('r [m]');
    ylabel('z [m]');
    axis equal
    axis([0 1 -1.2 1.2])
    set(gca, 'FontSize', 14);
    subplot(1, 3, 2);
    plot(VV(:, 1), VV(:, 2), '-k');
    hold on
    plot(SENSOR_TPRB.R, SENSOR_TPRB.Z, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4)
    title('MagProbe Position')
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
    subplot(1, 3, 3);
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

% % 再構成する磁束のサイズ
% % cm単位の範囲の指定
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

tic

% 逆問題を解いて求めたプラズマ電流およびコイルに流れる電流から、磁気面を計算
[PSI, DELGE, RCCS, ZCCS, XPSI] = INTER(PARAM, 1, GETA, CR, CZ, FFOUT, ExtCOIL, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, NINT);

toc

CCR = unique(CR);
CCZ = unique(CZ);
CCR(1) = [];
psi = reshape(PSI(1:numel(CCR) * numel(CCZ)), numel(CCZ), numel(CCR));

if CONFIG.ShowFig

    if CONFIG.DataType == 0
        v = linspace(-20, 20, 81);

        % figure
        % contour(CCR, CCZ, psi * 1000, v, '-k');
        % % contour(CCR, CCZ, (psi - imresize(REF.Flux, size(psi))) * 1000, v, '-k');
        % % image((psi - imresize(REF.Flux, size(psi))) * 1000, 'CDataMapping', 'scaled')
        % hold on
        % contour(REF.R, REF.Z, REF.Flux * 1000, v, '--m'); % 正解
        % contour(CCR, CCZ, psi, [0 0], 'LineColor', 'c', 'LineWidth', 2);
        % contour(REF.R, REF.Z, REF.Flux, [0 0], 'LineColor', 'm', 'LineWidth', 2);
        % plot(CCSDAT.RGI, CCSDAT.ZGI);
        % hold off
        % xlabel({'r (m)'});
        % ylabel({'z (m)'});
        % title("Reconstructed flux")
        % axis equal

        VV = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);
        figure
        hold on
        contour(CCR, CCZ, psi * 1000, v, '-k');
        contour(REF.R, REF.Z, REF.Flux * 1000, v, '--m'); % 正解
        contour(CCR, CCZ, psi, [0 0], 'LineColor', 'c', 'LineWidth', 2);
        contour(REF.R, REF.Z, REF.Flux, [0 0], 'LineColor', 'm', 'LineWidth', 2);
        scatter(SENSOR_TPRB.R, SENSOR_TPRB.Z);
        scatter(SENSOR_FLXLP.R, SENSOR_FLXLP.Z);
        plot(CCSDAT.RGI, CCSDAT.ZGI);
        plot(VV(:, 1), VV(:, 2), '-k'); % 容器壁
        hold off
        xlabel({'r (m)'});
        ylabel({'z (m)'});
        title("Reconstructed flux")
        axis equal
    elseif CONFIG.DataType == 1
        figure
        contour(CCR, CCZ, psi, '-k', 'LevelStep', 0.0003);
        hold on
        % contour(CCR, CCZ, psi, [0 0], 'LineColor', 'c', 'LineWidth', 2);
        plot(CCSDAT.RGI, CCSDAT.ZGI);
        hold off
        xlabel({'r (m)'});
        ylabel({'z (m)'});
        title("Reconstructed flux")
        axis equal
    end

    %% Data positions for 2D plot
    PLOT.R = 0.1:0.01:0.9;
    PLOT.Z = -1:0.02:1;
end

fclose('all');
% end

%% Main kokomade
