function MCCS_main_series_191030()
close all
clear all
fclose all;

% save figure?
SaveFig = 1

% 0: without, 1: average, 2: normalization between B and PSI,
% 3: normalization between B and PSI and between PSI_in and PSI_out
Normalization = 0;

% phase select (after merging: 1, 9950 us: 2, 9450 us: 3, 9530 us: 4, 8000 us: 5, 9510 us: 6, vacuum: 7)
select_phase = 6;

% L-curve iteration
Adopted_num = 50 : 1 : 55;
size_Adopted_num = size(Adopted_num);
assignin('base', 'Adopted_num', Adopted_num);

% input data
foldername_TF = '180528';
shotnum_TF = '001';
foldername = '180515';
shotnum = '010';

time_CCS = 9.525;
delTime = 0.005;
CCS_distance = 1000


f = fopen('CCS_center_position_180515010_modified.txt', 'r');
% f = fopen('CCS_center_position_180515009_modified.txt', 'r');
textscan(f,'%s',1,'Delimiter','\n');
for i = 1:31
    temp = textscan(f,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    time(i) = str2double(temp_m(1));
    CCSCenterZL(i) = str2double(temp_m(6));
    CCSCenterZU(i) = str2double(temp_m(7));
%     if 52 < i
%         CCSCenterZL(i) = CCSCenterZL(51);
%         CCSCenterZU(i) = CCSCenterZU(51);
%     end
end

figure
plot(time, CCSCenterZL)
hold on
plot(time, CCSCenterZU)

%%%%% Check here!!! %%%%%
% f1 = fopen('dummy.txt', 'w');
% f1 = fopen('ErrorBetweenReconstAnd2Darray_180515010_WithNormAndWeight_191120.txt', 'w');
% f1 = fopen('ErrorBetweenReconstAnd2Darray_180515010_WithAverageweight.txt', 'w');
% f1 = fopen('ErrorBetweenReconstAnd2Darray_180515010_WithNorm_191120.txt', 'w');
f1 = fopen('ErrorBetweenReconstAnd2Darray_180515010_Without_191120.txt', 'w');


fprintf(f1, 'time discrepancy_out discrepancy_in discrepancy_total AdptedNumForReconstruction\n');

    
TimeItrIdx_s = 1;
% TimeItrIdx_s = 49;
TimeItrIdx_e = 21;


hold on
plot(time(TimeItrIdx_s), CCSCenterZU(TimeItrIdx_s), 'o')
hold on
plot(time(TimeItrIdx_e), CCSCenterZU(TimeItrIdx_e), 'o')
% pause
for TimeItrIdx = TimeItrIdx_s : TimeItrIdx_e
    close all
            
    fprintf('iteration %i\n', TimeItrIdx)
    time_CCS = time(TimeItrIdx)

    xl = 0.35;
    xu = 0.35;

    [ITSKP,IT,AA,FF,NMAX,JMAX,n100,n50,NAPB,NFLX,NCCN,KNN,KSN,FC,FLXLP,BSNSR,GETA_YN,AUTO,AUTOINPUT,BZ_normalization_value,PSI_normalization_value] = CCS_UTST_L_curve_iteration_offline_data_input_190617_2(foldername_TF, shotnum_TF, foldername, shotnum, Adopted_num,xu,CCSCenterZU(TimeItrIdx),xl,CCSCenterZL(TimeItrIdx),time_CCS,select_phase,Normalization)

%     figure
%     subplot(1,2,1)
%     plot(FF(1:37))
%     subplot(1,2,2)
%     plot(FF(38:84))
%     pause
    
    for i = 1 : size_Adopted_num(2)
%         fclose('all');
        [C,W,U,V,FFOUT,XBFR,XMT,XGETA,GET] = SVD_MT_190623(ITSKP,IT,AA,FF,NMAX,JMAX,n100,n50,...
        0,0.0D0,NAPB,NFLX,NCCN,KNN,KSN,FC,FLXLP,BSNSR,GETA_YN,AUTO,AUTOINPUT, Adopted_num(i));
%     if AfterPhase == 0
%         CCS = 2;
%     else
%         CCS = 1;
%     end
        CCS = 2; % CCS�̐�
        NE = 3; % CCS�Z�N�V������
        size_FFOUT = size(FFOUT)

        half_norm = sqrt((sum((FFOUT).^2)));
        FF_reconst = AA(1 : 73 + CCS * NE * 3 -1, 1 : size_FFOUT(2)) * FFOUT';
        Error = FF_reconst - FF(1 : 73 + CCS * NE(1) * 3 -1)';
        Solution_norm = half_norm;
        Error_norm = sqrt(sum(Error.^2));

        SolNorm(i) = Solution_norm;
        ErrNorm(i) = Error_norm;
    end

    data_size = size(ErrNorm);
    for i = 2 : data_size(2) - 1 % �ŏ��ƍŌ�̓_�Ɋւ��Ă͋ȗ�����`�ł��Ȃ�
        vec1 = [ErrNorm(i - 1) - ErrNorm(i), SolNorm(i - 1) - SolNorm(i)]
        vec2 = [ErrNorm(i + 1) - ErrNorm(i), SolNorm(i + 1) - SolNorm(i)]
        kappa(i - 1) = acos(dot(vec1, vec2) / (norm(vec1) * norm(vec2))) * (180 / pi)
    end
    figure()
    plot(Adopted_num(2 : size_Adopted_num(2) - 1), kappa)
%     pause

    assignin('base','SolNorm',SolNorm);
    assignin('base','ErrNorm',ErrNorm);

    AdoptedNum = Adopted_num(2 : size_Adopted_num(2) - 1);
    [kappa_min, kappa_min_idx] = min(kappa);
    AdptedNumForReconstruction = AdoptedNum(kappa_min_idx);
    ANR(TimeItrIdx - TimeItrIdx_s + 1) = AdoptedNum(kappa_min_idx);

    [Error_2D_reconst_norm, CCR, CCZ, PSI_reconst, discrepancy, discrepancy_out, discrepancy_in] = MCCS_reconstruction_sub_2(select_phase, foldername, shotnum, foldername_TF, shotnum_TF, AdptedNumForReconstruction,xu,CCSCenterZU(TimeItrIdx),xl,CCSCenterZL(TimeItrIdx),time_CCS,Normalization,PSI_normalization_value)

    Discrepancy_out(TimeItrIdx - TimeItrIdx_s + 1) = discrepancy_out
    Discrepancy_in(TimeItrIdx - TimeItrIdx_s + 1) = discrepancy_in
    Discrepancy(TimeItrIdx - TimeItrIdx_s + 1) = discrepancy
    fprintf(f1, '%f %f %f %f %i\n', time_CCS, discrepancy_out, discrepancy_in, discrepancy, AdptedNumForReconstruction);

    assignin('base','PSI_reconst',PSI_reconst);
    assignin('base','CCR',CCR);
    assignin('base','CCZ',CCZ);

    f = figure()
    contourf(CCR , CCZ, PSI_reconst, 300);
    c = colorbar;
    c.Label.String = 'Poloidal flux [mWb]';
    caxis([-0.06 0.06]);
    xlabel({'R [m]'});
    ylabel({'Z [m]'});
    axis equal
    
    hold on
    MASK1 = dlmread('MASK1.txt');
    fill(MASK1(:,1),MASK1(:,2),'w')
    hold on
    fill(MASK1(:,1),-MASK1(:,2),'w')
    axis equal
    VV = dlmread('VacuumVesselMeshPoints.txt');
    PLOTW = 1
    if (PLOTW == 0)
        axis([min(VV(:,1)) max(VV(:,1)) -0.5 0.5])%min(VV(:,2)) max(VV(:,2))])
    elseif (PLOTW == 1)
        axis([min(VV(:,1)) max(VV(:,1)) min(VV(:,2)) max(VV(:,2))])%min(VV(:,2)) max(VV(:,2))])
        view(2)
    end
    
    
%     save('PSI_reconst_am_1911010', PSI_reconst)
    
%     pause

    if SaveFig == 0
        FigName = strcat('FIGS_191120/Flux_surface_', num2str(round(1000 * time_CCS)))
        FigName = strcat(FigName, '.png')
        saveas(f, FigName);
    end

%     time_CCS = time_CCS + delTime;
end % if TimeItrIdx == 1

figure()
plot(Discrepancy)
pause

figure()
plot(ErrNorm, SolNorm, 'o')
pause

figure()
for i = 1 : size_Adopted_num(2)
    plot(ErrNorm(i), SolNorm(i), 'o', 'DisplayName', num2str(Adopted_num(i)))
    hold on
end
legend










