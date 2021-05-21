function [Error_2D_reconst_norm, CCR, CCZ, PSI_reconst, discrepancy, discrepancy_out, discrepancy_in] = MCCS_reconstruction_sub_2(select_phase, foldername, shotnum, foldername_TF, shotnum_TF, Adopted_num,xu,yu,xl,yl,time_CCS,Normalization,PSI_normalization_value)
    % *****************************************************************
    % *  CAUCHY CONDITION SURFACE METHOD  (quadratic element version) *
    % *****************************************************************
%     close all
%     fclose('all');
    format long e

% Adopted_num = 44;    
    
% ������͂��Ă�����
% foldername_TF = '180528';
% shotnum_TF = '001';
% foldername = '180515';
% shotnum = '010';

% EF���ӂ��Ď�����������
% TF�V���b�g�d����1.5�A�ʏ���d��2.7kV�ŏ[�d���Ă��܂������߁A190128001��TF�V���b�g�i2.7kV�j���Q��
% ��002,003,004��TF��0kV�i���̃R�[�h���ŉ񂷂Ȃ�TF�����d���Ȃ���΂Ȃ�Ȃ��j
% foldername_TF = '190128';
% shotnum_TF = '001';
% foldername = '190125';
% shotnum = '007';

% �^����d��������
% foldername_TF = '190128';
% shotnum_TF = '001';
% foldername = '190128';
% shotnum = '002';

EF_voltage = 120;
        
    % ���̂�����̎��O���蓖�Ăɂ�閳�ʂȔz��̍쐬���v�Z���Ԃ�x�����Ă���
    nsemx = 2052;   % Max. No. of sensors (FLAXLOOP:1024,T-Probe:1024)
    n100 = 250;   % Max. No. of given data ( > NAPB+NFLX+MXCCS)%250
    n50 = 150;     % Max. No. of unknowns ( > 2*MXCCS)    
    Nedp = 150;
    AA = zeros(n100,n50); 
    FF = zeros(1,n100);
    FFOUT = zeros(1,n50);  
    GHR = zeros(1,300); 
    GHZ = zeros(1,300); 
%
    ITEND = 0; % ushiki
    EPS1 = 0; % ushiki
    ECIGRP = 0;% ushiki
    I1I0 = 0;% ushiki
    PSI = 0;% ushiki
% 
    fid222 = fopen('@extendedGII.txt','w');
    myformat = '%d %d\r\n';
    X1 = -1.0;
    X2 = 0.0;
    X3 = 1.0;
    IMAX = 20;
    CMAX = IMAX;
    GEND = 1000.0;
    GSTAT = -GEND;
    DEL = (GEND-GSTAT)/CMAX;
for I = 1:IMAX+1
    CIM1 = I-1;
    GII = GSTAT + CIM1*DEL;
    F1=GII*(GII-1.0)*0.5;
    F2=1.0-GII^2;
    F3=GII*(GII+1.0)*0.5;
    XCO=X1*F1+X2*F2+X3*F3;
    fprintf(fid222,myformat,GII,XCO);
end

fprintf('wahahahahahaha!!!!\n'); % IPR
WAHAHA = fopen('WAHAHA.txt','w');% IPR
RMYU0=4.0*pi*1.0D-7;
%CHG = pi/180;
MXCCS = 20;          % MAX! / NUMBER OF ELEMENTS ON THE CCS
MXINT = 10000;       % MAX! / NUMBER OF INTERNAL POINTS
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
% prompt = 'automatic or manual �iCCS-2�j? = (1/0)\n';
% AUTO = input(prompt);
AUTO = 1;
% if (AUTO == 1)
    prompt = 'merging / after merging / merging(9550us) / merging(9450us) / merging(9530us) / vacuum(8000us) / merging(9510us) / vacuum shot = (0/1/2/3/4/5/6/7)\n';
%     select_phase = input(prompt);
%     select_phase = 6;
    if (select_phase == 0)
        AUTOINPUT = load('INPUT/INPUT_merging.dat');
    elseif (select_phase == 1)
        AUTOINPUT = load('INPUT/INPUT_after_merging.dat');
        time_CCS = 9.65;
    elseif (select_phase == 2)
        AUTOINPUT = load('INPUT/INPUT_merging_9550.dat');
        time_CCS = 9.55;
    elseif (select_phase == 3)
        AUTOINPUT = load('INPUT/INPUT_merging_9450.dat');
        time_CCS = 9.45;
    elseif (select_phase == 4)
        AUTOINPUT = load('INPUT/INPUT_merging_9530.dat');
%         time_CCS = 9.53;
    elseif (select_phase == 5)
        AUTOINPUT = load('INPUT/INPUT_vacuum_8000.dat');
        time_CCS = 8.00;
    elseif (select_phase == 6)
        AUTOINPUT = load('INPUT/INPUT_merging_9510.txt');
%         time_CCS = 9.51;
    elseif (select_phase == 7)
        AUTOINPUT = load('INPUT/INPUT_vacuum_shot.dat');
        prompt = 'input the time [ms]\n';
        time_CCS = input(prompt);
    elseif (select_phase == 8)
        AUTOINPUT = load('INPUT/INPUT_Xie.txt');
        time_CCS = time_input;
    end
% else
%     AUTOINPUT = 0;
% end
% if (AUTO == 0)
%     prompt = 'JT-60U? / UTST / UTST merging / UTST breakdown / UTST merging 2 / UTST merging unequal  ? =(0/1,2,3,4,5/6,7/8/9/10,..)\n';
%     IUTST = input(prompt);
% else
IUTST = AUTOINPUT(1);
% end

eddy_time = 0;
SENSio = 0;
LWALL = 0;
% if (IUTST == 8)
%     prompt = 'Choose time (7.1/7.2/7.3/7.4/7.5/7.6/7.7/7.8msec)\n';    
%     eddy_time =input(prompt)*1000;
%     prompt = 'Using sensor only inside? / Using sensor outside? =(0/1)\n';        
%     SENSio =input(prompt);
%     prompt = '�C���{�[�h���̃Z���T�[�g��? / �g��Ȃ�? =(0/1)\n';        
%     LWALL =input(prompt);
% end

[CXK,CYK,SOLK,MSK] = FNGRPH_UTST(IUTST,eddy_time);
fprintf('*** Start of PreUTST\n');
[RS,ZS,ITYPE,TET,select_phase, Plasma_Current_minus_offset, Psi_z_nega_200, Psi_z_nega_100, Psi_z0, Psi_z100, Psi_z200, R_2d] = PreUTST(IUTST,WAHAHA,nsemx,eddy_time,SENSio,LWALL,select_phase, time_CCS, foldername_TF, shotnum_TF, foldername, shotnum, EF_voltage); % OK
% plot sensor position
VV = dlmread('VacuumVesselMeshPoints.txt');
SEN0 = dlmread('SENPOS0.txt');
SEN1 = dlmread('SENPOS1.txt');
% figure('Name','Sensor Position','NumberTitle','off')
% plot(VV(:,1),VV(:,2),'-k', SEN0(:,1),SEN0(:,2),'o','MarkerEdgeColor','b',...
%     'MarkerFaceColor','b','MarkerSize',5)
% hold on
% plot(SEN1(:,1),SEN1(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r',...
%     'MarkerSize',4)
% xlabel({'r [m]'});
% ylabel({'z [m]'});
% axis equal
% fprintf('*** End of PreUTST\n');   
% axis([0 1 -1.2 1.2])
% set(gca, 'FontSize',14);
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX   

% 
if (AUTO == 0)
    prompt = 'Use Geta iteration? (Yes/No)=(1/0)\n';
    GETA_YN = input(prompt);
%
    prompt = 'Conforming EddyCrnt BE? / Non-Conforming BE? =(0/1)\n';
    NONC = input(prompt);
else
    GETA_YN = AUTOINPUT(2); %���ʔ������邩�ǂ���
    NONC = AUTOINPUT(3);
end
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
[CR,CZ,FFDAT,R0,Z0,RR,NCOUNT,RC,ZC,ECI,NTPB,NNPB,NAPB,NFLX,FLXLP,BSNSR,CAPPER,...
REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,IDECCS,SOU,REVN,ZEVN,GETA,NE,KCMX,ISHT,NINT,...
R_plot,Z_plot,CCS,MSEC] = INPUT(n100,Nedp,IUTST,MXINT,WAHAHA,NONC,AUTO,AUTOINPUT,xu,yu,xl,yl); % OK


%
% plot CV segment
COIL = dlmread('@UTST_CoilGeom.txt');
VVMESH = dlmread('VacuumVesselMeshPoints.txt');
VVNODE = dlmread('VacuumVesselNodePoints.txt');
VVSEG = dlmread('VacuumVesselSegments.txt');

figure('Name','CVsegment','NumberTitle','off')
plot(VV(:,1),VV(:,2),'-k')
hold on
plot(VVMESH(:,1),VVMESH(:,2),'bx','MarkerSize',4)
hold on
plot(VVNODE(:,1),VVNODE(:,2),'ko','MarkerFaceColor','k','MarkerSize',3)
hold on
plot(COIL(:,1), COIL(:,2), 's', 'MarkerSize', 3)
% hold on1
% plot(RES(1:3),ZES(1:3),'ko',RES(1:3),ZES(1:3),'-k')
xlabel({'r [m]'});
ylabel({'z [m]'});
axis equal
fprintf('*** End of PreUTST\n');   
axis([0 1 -1.2 1.2])
set(gca, 'FontSize',14);
% pause

fprintf(WAHAHA,'Check of FFDAT(I)\n');   
%for I = 1:NAPB+NFLX
fprintf(WAHAHA,'%d %d\r\n',horzcat((1:NAPB+NFLX)',FFDAT(1:NAPB+NFLX)'));   
%end
%
TOTMOD = 1.0; 
NCCS(1:CCS) = NE(1:CCS)*2;
%
ICONT = 0;
%  9999 CONTINUE ! ICONT�̍X�V�iCCS����蒼���j!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
if (ICONT == 0)  % ICONT>0 �Ȃ�CCS����蒼���B
    if (IDECCS > 0) %! �ȉ~�^CCS or D�^CCS
       %%
        %*******************************************************
        %  'D'�^CCS                                         ****     
        %*******************************************************
        [RCCS,ZCCS,RCCN,ZCCN,NCCN,NE,NCCS] = D_CCS(R0,Z0,RR,CAPPER,MXCCS,CCS,MSEC);
    else
        RCCS = zeros(CCS,MXCCS+1);
        ZCCS = zeros(CCS,MXCCS+1);%1
        RCCN = zeros(CCS,MXCCS+1);
        ZCCN = zeros(CCS,MXCCS+1);
        NCCN = zeros(1,CCS);%ushiki
        % C*******************************************************
        % C  �ȉ~�^CCS                                        ****     
        % C*******************************************************
        for III = 1:CCS
            DTHETA = 2.0*pi/NCCS(III);
            TRIG=0.0;      %�O�p�x  
%        fprintf(WAHAHA,'%d %d\r\n',I,FFDAT(I));   
        %        WRITE(IPR,3)
        %     3 FORMAT(1X,'Meshpoints for CCS',/,
        %      &       4X,'No.',9X,'R',14X,'Z') 
            for I = 1:NCCS(III)
                THETA = pi/2.0-DTHETA*(I-1);  % CCS�ł̐ϕ��́A���v���
                RCCS(III,I) = R0(III) + RR(III)*cos(THETA+asin(TRIG)*sin(THETA));
                ZCCS(III,I) = Z0(III) + CAPPER(III)*RR(III)*sin(THETA);
                fprintf(WAHAHA,'%d %d %d\r\n',I,RCCS(III,I),ZCCS(III,I));   
            end
        end
    end
else
    % *******************************************************
    %    LCMS�ɑ�����CCS    (ICONT>0 �̂Ƃ�)              ****
    % *******************************************************
    %       WRITE(IPR,3)
    for III = 1:CCS
        fprintf(WAHAHA,'%s\r\n',' ');
        for I = 1:NCCS(III)
            II=NCCS(III)+2-I;
            if (II > NCCS(III))
                II=1;
            end
            RCCS(III,I) = GHR(II);
            ZCCS(III,I) = GHZ(II);
            fprintf(WAHAHA,'%d %d %d\r\n',I,RCCS(III,I),ZCCS(III,I));   
        end 
    end
end  
% 
% -----------------------------------------------------------------------
%  CCS��̔�K���v�f�ߓ_���W�̍쐬
% 
fid12 = fopen('MeshPoints.txt','a'); 
fid13 = fopen('DiscontinuousNodePoints.txt','w');
fid14 = fopen('CurvCCS.txt','w');
frewind(fid14);
fid15 = fopen('CurvCCS_Final.txt','w'); 
frewind(fid15);
myformat = '%d %d\r\n';
for III = 1:CCS
    RCCS(III,NCCS+1)=RCCS(III,1);
    ZCCS(III,NCCS+1)=ZCCS(III,1);
    NCCN(III)=NE(III)*3;
    I = 1:NE(III);
    RCCN(III,3*I-2) = (5.*RCCS(III,2*I-1) + 5.*RCCS(III,2*I)-RCCS(III,2*I+1))/9;
    ZCCN(III,3*I-2) = (5.*ZCCS(III,2*I-1)+5.*ZCCS(III,2*I)-ZCCS(III,2*I+1))/9;
    RCCN(III,3*I-1) = RCCS(III,2*I);
    ZCCN(III,3*I-1) = ZCCS(III,2*I);
    RCCN(III,3*I)   = (5.*RCCS(III,2*I+1)+5.*RCCS(III,2*I)-RCCS(III,2*I-1))/9;
    ZCCN(III,3*I)   = (5.*ZCCS(III,2*I+1)+5.*ZCCS(III,2*I)-ZCCS(III,2*I-1))/9;
    %
    for I = 1:NCCS(III)
        fprintf(fid12,myformat,RCCS(III,I),ZCCS(III,I));
    end
    fprintf(fid12,myformat,RCCS(III,1),ZCCS(III,1));
    %fprintf(fid12,'%s\r\n','= =');
    fprintf(WAHAHA,'%s\r\n',' Node points for CCS');
    for I = 1:NCCN(III)
        fprintf(fid13,myformat,RCCN(III,I),ZCCN(III,I)); 
        fprintf(WAHAHA,myformat,RCCN(III,I),ZCCN(III,I));
    end
%    fprintf(fid13,'%s\r\n','= ='); 
% 
    %  Draw a fine curve to express the shape of CCS 
    MAXG = 200;  
    GMAX = MAXG;
    DEL = 2.0/GMAX;
    for I = 1:NE(III)
        for J = 1:MAXG+1
            CJM1 = J-1;
            GII = -1.0+DEL*CJM1;
            F1 = GII*(GII-1.0D0)*0.5;
            F2 = 1.0D0-GII^2;
            F3 = GII*(GII+1.0D0)*0.5;
            RGI = RCCS(III,2*I-1)*F1+RCCS(III,2*I)*F2+RCCS(III,2*I+1)*F3;
            ZGI = ZCCS(III,2*I-1)*F1+ZCCS(III,2*I)*F2+ZCCS(III,2*I+1)*F3;
            fprintf(fid14,myformat,RGI,ZGI); 
            fprintf(fid15,myformat,RGI,ZGI); 
        end
    end
% CCSMESH = dlmread('@D_CCS_CheckWrite.txt');
% CCSNODE = dlmread('DiscontinuousNodePoints.txt');
% hold on
% plot(CCSMESH(:,1),CCSMESH(:,2),'-k', CCSMESH(:,1),CCSMESH(:,2),'bx',CCSNODE(:,1),CCSNODE(:,2),'ko',...
%      VVSEG(:,1),VVSEG(:,2),'r+','MarkerSize',2)
%    fprintf(fid14,'%s\r\n','= ='); 
% 
RCCSSP(III,:) = spline(1:NCCN+1,horzcat(RCCN(III,1:NCCN),RCCN(III,1)),1:1/10:NCCN+1);
ZCCSSP(III,:) = spline(1:NCCN+1,horzcat(ZCCN(III,1:NCCN),ZCCN(III,1)),1:1/10:NCCN+1);
hold on
plot(RCCSSP(III,:),ZCCSSP(III,:),'-k')
hold on
plot(RCCS(III,1:NCCS),ZCCS(III,1:NCCS),'bx','MarkerSize',4)
hold on
plot(RCCN(III,1:NCCN),ZCCN(III,1:NCCN),'ko','MarkerFaceColor','k','MarkerSize',3)
end
% �Z���T�M��FFDAT��FF�ɏ���������
FF(1:NAPB+NFLX) = FFDAT(1:NAPB+NFLX);
FF(NAPB+NFLX+1:NAPB+NFLX+sum(NCCN))=0.0;
% prompt = 'Constraint of total plasma current? =(0/1)\n'; %�v���Y�}�d������
% ipconst = input(prompt);
ipconst = 0;

if (ipconst == 1)
%        FF(NAPB+NFLX+sum(NCCN)+1)=-66410*RMYU0;
    % Ip���S�X�L�[�̃f�[�^���C���v�b�g
    FF(NAPB+NFLX+sum(NCCN)+1) = Plasma_Current_minus_offset(round(time_CCS / (0.5 * 0.001))) * 1000 * RMYU0;
end
% 
%  ******************************************************************************************************
% 
%      iteration Start
% 
%CC ITMX=1;
ITMX = 1000;
DELGE = 0.0;
OLDEL = 1.0;
GETA = 0.0; % ! ���ʂ̏����l���[���ɂ���
OLDGT = GETA;
fid90 = fopen('GETA.txt','w'); %90
OMGA = 1.9;
%
% ***************************
% ***************************
if (AUTO == 0)
    fprintf('���ʃT�[�`�� INTER�̔����ł����H (0)\n');
    prompt ='���ʃT�[�`��SVD_MT�����݂̂ł��H (1)\n';
    ITSKP = input(prompt);
else
    ITSKP = AUTOINPUT(27);
end
if (ITSKP > 0)
    ITMX=1;
fprintf('���ʃT�[�`�� (1) SVD_MT�����݂̂ł��܂�\n');
else
fprintf('���ʃT�[�`�� (0) INTER�̔����ł����܂�\n');
end
% ***************************
% ***************************
%       
for IT = 1:ITMX
    GETA = GETA + DELGE;
    GETA = OMGA*GETA+(1.0D0-OMGA)*OLDGT;
    fprintf('//\n');
    fprintf('//\n');
    fprintf('****************************************************\n');
    fprintf('%s%d%s\r\n','** (IT) Iteration count =',IT,' ********************');
    fprintf('%s %d  %d\r\n','** (GETA/DELGE) =',GETA,DELGE);
    %******************************************************************************************************
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    [FC,BR,BZ,PSIFLX,PSIC,AMYU0,AA,FF,GETA,NCCS,fid99,fid100] = FORM(WAHAHA,IT,OLDGT,...
    GETA,FF,AA,n100,n50,ECI,KCMX,RS,ZS,RC,ZC,ITYPE,TET,NTPB,NNPB,NAPB,NFLX,...
    RCCN,ZCCN,NCCN,NCCS,RCCS,ZCCS,REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,Nedp,NONC,MXCCS,...
    RMYU0,NE,CCS,ipconst);%NE uhiski OK

    %
    % Solve Matrix Equation
    % Singular Value Decompoition
    if (ipconst == 1)
        NMAX = NAPB + NFLX + sum(NCCN)+1;% ushikiip
    else
        NMAX = NAPB + NFLX + sum(NCCN);% ushikiip
    end
    JMAX = sum(NCCN) + sum(NCCN) + sum(KNN) + sum(KSN);
    fprintf('NCCN/KNN/KSN = %d %d %d\r\n',sum(NCCN),sum(KNN),sum(KSN));
    
    % �d�ݕt���ŏ����@�������i�ߓ��j
    % AA�F�W���s��AFF�F����l�x�N�g���i���C�Z���T39�A�����Z���T33�j
    % AA*FFOUT2=FF
    size(AA(1:81, 1:60))
    Coeff_mat = AA(1:73 + CCS * NE * 3 -1, 1:60);
    % �R�C���d����^����������
    Sig_vec = transpose(FF(1:73 + CCS * NE(1) * 3 -1))
    % �d�݃x�N�g���͐��łȂ���΂Ȃ�Ȃ�
    Sig_abs = ((Sig_vec.^2).^(0.5))
    MF_abs_ave = sum(Sig_abs(1 : 39)) ./ 39
    FX_abs_ave_in = sum(Sig_abs(40 : 58)) ./ 18
    FX_abs_ave_out = sum(Sig_abs(59 : 72)) ./ 16
    FX_CCS_abs_ave = sum(Sig_abs(73 : 73 + CCS * NE(1) * 3 -1)) ./ (CCS * NE(1) * 3)
    
    NV_b = max(abs(FF(1 : 39)))
    NV_f = max(abs(FF(40 : 73 + CCS * NE(1) * 3 -1)))
    PSI_normalization_value = PSI_normalization_value / max(PSI_normalization_value);
    PSI_normalization_value_in = mean(PSI_normalization_value(1, 1 : 19))
    PSI_normalization_value_out = mean(PSI_normalization_value(1, 20 : 35))
    PSI_normalization_value_ratio = (PSI_normalization_value_out / PSI_normalization_value_in)^(0.5)
    
    if Normalization == 3
        FF(1 : 39) = FF(1 : 39) / NV_b;
        FF(40 : 73 + CCS * NE(1) * 3 -1) = FF(40 : 73 + CCS * NE(1) * 3 -1) / NV_f;
        for i = 1 : 39
            AA(i, :) = AA(i, :) / NV_b;
            plot(AA(i, :))
        end
        for j = 40 : (73 + CCS * NE(1) * 3 -1)
            if and(39 < j, j < 58)
                AA(j, :) = AA(j, :) / NV_f * PSI_normalization_value_ratio;
                FF(j) = FF(j) * PSI_normalization_value_ratio;
            elseif and(57 < j, j < 73)
                AA(j, :) = AA(j, :) / NV_f;
                FF(j) = FF(j);
            elseif and(72 < j, 73 + CCS * NE(1) * 3 -1)
                AA(j, :) = AA(j, :) / NV_f;
                FF(j) = FF(j);
            end
        end
    elseif Normalization == 2
        FF(1 : 39) = FF(1 : 39) / NV_b;
        FF(40 : 73 + CCS * NE(1) * 3 -1) = FF(40 : 73 + CCS * NE(1) * 3 -1) / NV_f;
        for i = 1 : 39
            AA(i, :) = AA(i, :) / NV_b;
            plot(AA(i, :))
        end
        for j = 40 : (73 + CCS * NE(1) * 3 -1)
            if and(39 < j, j < 58)
                AA(j, :) = AA(j, :) / NV_f;
                FF(j) = FF(j);
            elseif and(57 < j, j < 73)
                AA(j, :) = AA(j, :) / NV_f;
                FF(j) = FF(j);
            elseif and(72 < j, 73 + CCS * NE(1) * 3 -1)
                AA(j, :) = AA(j, :) / NV_f;
                FF(j) = FF(j);
            end
        end
    elseif Normalization == 1
        AA(1 : 26, :) = AA(1 : 26, :) ./ (MF_abs_ave)
        AA(27 : 31, :) = AA(27 : 31, :) ./ (MF_abs_ave)
        AA(32 : 39, :) = AA(32 : 39, :) ./ (MF_abs_ave)
        AA(40 : 57, :) = AA(40 : 57, :) ./ (FX_abs_ave_in)
        AA(58 : 72, :) = AA(58 : 72, :) ./ (FX_abs_ave_out)
        AA(73 : 73 + CCS * NE(1) * 3 -1, :) = AA(73 : 73 + CCS * NE(1) * 3 -1, :) ./ (FX_CCS_abs_ave)
        FF(1 : 26) = FF(1 : 26) ./ (MF_abs_ave)
        FF(27 : 31) = FF(27 : 31) ./ (MF_abs_ave)
        FF(32 : 39) = FF(32 : 39) ./ (MF_abs_ave)
        FF(40 : 57) = FF(40 : 57) ./ (FX_abs_ave_in)
        FF(58 : 72) = FF(58 : 72) ./ (FX_abs_ave_out)
        FF(73 : 73 + CCS * NE(1) * 3 -1) = FF(73 : 73 + CCS * NE(1) * 3 -1) ./ (FX_CCS_abs_ave)
    else
        AA = AA;
        FF = FF;
    end
   
    %  Modified Truncated SVD for Smoothness
%     Adopted_num = 43;
    [C,W,U,V,FFOUT,XBFR,XMT,XGETA,GET] = SVD_MT_190502(ITSKP,IT,AA,FF,NMAX,JMAX,n100,n50,...
    0,0.0D0,NAPB,NFLX,NCCN,KNN,KSN,FC,FLXLP,BSNSR,GETA_YN,AUTO,AUTOINPUT, Adopted_num);
    %
    % 
    fprintf(WAHAHA,'%s\r\n','  F  VECTOR at the end of iteration');
    fprintf(WAHAHA,'%s\r\n','Results along CCS');
    for N = 1:sum(NCCN)
        fprintf(WAHAHA,'N=   %d d��/dn=   %d',N,FFOUT(N));
    end
    for N = 1+sum(NCCN):2*sum(NCCN)
        fprintf(WAHAHA,'N=   %d    �� =   %d',N-sum(NCCN),FFOUT(N));
    end
%���m����
half_norm = sqrt((sum((FFOUT).^2)));
fprintf('%s%d\r\n','norm of the solution vector = ',half_norm);

size_FFOUT = size(FFOUT)
FF_reconst = AA(1 : 73 + CCS * NE * 3 -1, 1 : size_FFOUT(2)) * FFOUT';

% figure('Name','Check reconstructed value','NumberTitle','off')
% plot(FF, 'ko')
% hold on
% plot(FF_reconst, 'ro')

Error = FF_reconst - FF(1 : 73 + CCS * NE(1) * 3 -1)'
% figure('Name','Error','NumberTitle','off')
% plot(Error)

fprintf('Ax = b\nNorm of solution vector |x|')
half_norm
fprintf('Norm of error vector |Ax - b|')
% size(FF_reconst)
sqrt(sum(Error.^2))


% fprintf('Checking the accuracy of the reconstruction')
% fprintf('Between reconstructed and measurement value')
% size(FFOUT)
% size(Coeff_mat)
% size(Sig_vec)
% Reconst_at_sensor = Coeff_mat * FFOUT'
% figure('Name','reconstructed and measurement value','NumberTitle','off')
% plot(Reconst_at_sensor)
% hold on1
% plot(Sig_vec)


% 
% ***************************
% ***************************
    if (ITSKP > 0)
        GETA = XGETA;
        break
    else
    end
% C***************************
% C***************************

    [PSI,DELGE,RCCS,ZCCS,XPSI] = INTER(0,GETA,RCCS,ZCCS,CR,CZ,FFOUT,n100,n50,RC,ZC,...
    ECIGRP,ECI,KCMX,NAPB,NFLX,FLXLP,BSNSR,RS,ZS,TET,RCCN,ZCCN,...
    NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,WAHAHA,...
    MXCCS,MXINT,NCCS,NINT,RMYU0,NE,CCS);
% 
    EPSDEL = abs((DELGE-OLDEL)/OLDEL);
    if and(IT > 2, EPSDEL < 0.1)
        OMGA = OMGA + 5.0;
        fprintf('%s %d\r\n','    ******  OMGA changed to', OMGA);
    end
    OLDEL = DELGE;
    OLDGT = GETA;
    EPS1 = abs(DELGE/GETA);
    ITEND = IT;
    if (EPS1 < 1.0D-5)
        break
    end
end	
% 
%  ******************************************************************************************************
%      iteration End
%  ******************************************************************************************************
fprintf(' \n');

fprintf(' \n');
fprintf('********* Congraturation!!  GETA iteration Converged.\r\n');
fprintf('%s%d\r\n','****** No. of iterations required = ',ITEND);
fprintf('%s%d %d\r\n','** (GETA/DELGE) = ',GETA,DELGE);
fprintf('%s%d\r\n','** EPS=DABS(DELGE/GETA)= ',EPS1);
% 
[n50] = EDDYP(FFOUT,n50,KNE,KNN,NCCN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,REVN,ZEVN);
% plot eddy current
DISF = dlmread('EddyCurrentProfile.txt');
figure('Name','Eddy Current Plofile','NumberTitle','off')
plot(DISF(:,1),DISF(:,2),'-ko','MarkerEdgeColor','r','MarkerFaceColor','r',...
                                                             'MarkerSize',2) %�Q�d���č\��
xlabel({'Distance (m)'});
ylabel({'Eddy Current Density (MA/m)'});
axis([0 DISF(end,1) -0.2 0.2])
set(gca, 'FontSize',14);
dlmwrite(strcat('UTST_',num2str(IUTST),'/jeddyrec.txt'),DISF);
if and(IUTST >= 4,IUTST <= 5)
    axis([0 DISF(end,1) -0.2 0.2])
    set(gca, 'FontSize',14);
    hold on
    JEDDY = dlmread(strcat('UTST_',num2str(IUTST),'/jeddy3b.txt'));    
    plot(JEDDY(:,1),JEDDY(:,2),'-b')
    hold on
    JEDDY = dlmread(strcat('UTST_',num2str(IUTST),'/jeddy3c.txt'));    
    plot(JEDDY(:,1),JEDDY(:,2),'-r')
elseif (IUTST == 6)
    axis([0 DISF(end,1) -0.2 0.2])
    set(gca, 'FontSize',14);
    hold on
    JEDDY = dlmread(strcat('UTST_',num2str(IUTST),'/jeddy4c.txt'));    
    plot(JEDDY(:,1),JEDDY(:,2),'-b')
elseif (IUTST == 7)
    axis([0 DISF(end,1) -0.2 1.5])
    set(gca, 'FontSize',14);
    hold on
    JEDDY = dlmread(strcat('UTST_',num2str(IUTST),'/jeddy5c.txt'));    
    plot(JEDDY(:,1),JEDDY(:,2),'-b')
elseif (IUTST == 8)
    axis([0 DISF(end,1) -0.4 0.4])
    set(gca, 'FontSize',14);
    hold on
    JEDDY = dlmread(strcat('UTST_',num2str(IUTST),'/jeddy',num2str(IUTST),'_',num2str(eddy_time),'.txt'));    
    plot(JEDDY(:,1),JEDDY(:,2),'-b')
elseif or(IUTST == 9,IUTST == 10)
    axis([0 DISF(end,1) -0.4 0.4])
    set(gca, 'FontSize',14);
    hold on
    JEDDY = dlmread(strcat('UTST_',num2str(IUTST),'/jeddy',num2str(IUTST),'_new2.txt'));  
%     plot(JEDDY(:,1),JEDDY(:,2),'-g')
end
if (AUTO == 0)
    fprintf(' \n');
    fprintf('Time consuming, but you draw a whole flux profile?\n');
    prompt ='(Yes/No)=(1/0)\n';
    LXMAP = input(prompt);
else
    LXMAP = AUTOINPUT(33);
end
if (LXMAP > 0)
% �����̓��������v���b�g
[PSI,DELGE,RCCS,ZCCS,XPSI] = INTER(1,GETA,RCCS,ZCCS,CR,CZ,FFOUT,n100,n50,RC,ZC,...
ECIGRP,ECI,KCMX,NAPB,NFLX,FLXLP,BSNSR,RS,ZS,TET,RCCN,ZCCN,...
NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,WAHAHA,...
MXCCS,MXINT,NCCS,NINT,RMYU0,NE,CCS);
% Comparison_FieldSignal
COMPf = dlmread('Comparison_FieldSignal.txt');
COMPf2 = dlmread('Comparison_TotField.txt');
COMPfl = dlmread('Comparison_FluxSignal.txt');
COMPfl2 = dlmread('Comparison_TotFlux.txt');

% figure('Name','Check Reconstructed Field Value','NumberTitle','off')
% plot(COMPf2(:, 1), 'ko')
% hold on
% plot(COMPf2(:, 2), 'ro')

Error_field = (COMPf2(:,1)-COMPf2(:,2)).^2;
% figure('Name','Error (field)','NumberTitle','off')
% plot(Error_field)

% figure('Name','Check Reconstructed Flux Value','NumberTitle','off')
% plot(COMPfl2(:, 1), 'ko')
% hold on
% plot(COMPfl2(:, 2), 'ro')

Error_flux = (COMPfl2(:,1)-COMPfl2(:,2)).^2;
% figure('Name','Error (flux)','NumberTitle','off')
% plot(Error_flux(1:18))

fprintf('flux error (inboard) |g - g*| \n%d', sum(Error_flux(1:18)) / 18)
fprintf('\n flux error (all) |g - g*| \n%d', sum(Error_flux) / 33)
fprintf('\n field error |g - g*| \n%d', sum(Error_field) / 39)
fprintf('\n total error |g - g*| \n%d', (sum(Error_field) + sum(Error_flux)) / 72)

Flux_Reconst_at_sensor = COMPfl2(:,2);
Flux_Measured_at_sensor = COMPfl2(:,1);
Flux_size = size(Flux_Reconst_at_sensor);
Normalized_Flux_Error = zeros(Flux_size(1), 1);
for i = 1 : Flux_size
    Normalized_Flux_Error(i, 1) = abs((Flux_Reconst_at_sensor(i) - Flux_Measured_at_sensor(i)) / Flux_Measured_at_sensor(i));
end


% figure('Name','Normalized Flux Error','NumberTitle','off')
% plot(Flux_Measured_at_sensor, Normalized_Flux_Error, 'ko')


Field_Reconst_at_sensor = COMPf2(:,2);
Field_Measured_at_sensor = COMPf2(:,1);
Field_size = size(Field_Reconst_at_sensor);
Normalized_Field_Error = zeros(Field_size(1), 1);
for i = 1 : Field_size
    Normalized_Field_Error(i, 1) = abs((Field_Reconst_at_sensor(i) - Field_Measured_at_sensor(i)) / Field_Measured_at_sensor(i));
end


% figure('Name','Normalized Field Error','NumberTitle','off')
% plot(Field_Measured_at_sensor, Normalized_Field_Error, 'ko')


figure('Name','Comparison Field & Flux Signal','NumberTitle','off')
deg45 = dlmread('Line45deg.txt');
field_error = sqrt((sum((COMPf2(:,1)-COMPf2(:,2)).^2)) /numel(COMPf2(:,1)));
fprintf('%s%d\r\n','error of the field sensors = ',field_error);
% subplot(1,2,1);plot(COMPf(:,1),COMPf(:,2),'ko',COMPf2(:,1),COMPf2(:,2),'bo',deg45(:,1),deg45(:,2),'-r')
% xlabel({'Given field signal (T)'});
% ylabel({'Reconstructed field signal (T)'});
% axis equal
% axis([-0.25 0.25 -0.25 0.25])
% set(gca, 'FontSize',14);
%
% Comparison_FluxSignal
flux_error = sqrt((sum((COMPfl2(:,1)-COMPfl2(:,2)).^2)) /numel(COMPfl(:,1)));
fprintf('%s%d\r\n','error of the flux sensors = ',flux_error);
% subplot(1,2,2);plot(COMPfl(:,1),COMPfl(:,2),'ko',COMPfl2(:,1),COMPfl2(:,2),'bo',deg45(:,1),deg45(:,2),'-r')
% xlabel({'Given flux signal (Wb/rad)'});
% ylabel({'Reconstructed flux signal (Wb/rad)'});
% axis equal
% % axis([-0.005 0.025 -0.005 0.025])
% axis([-0.05 0.05 -0.05 0.05])
% set(gca, 'FontSize',14);
%
%



if (AUTO == 0)
    prompt = 'Plot a whole or zoomed flux profile? (1/0)\n';
    PLOTW = input(prompt);
else
    PLOTW = AUTOINPUT(34);
end
%         figure('Name','Flux Profile','NumberTitle','off')
%        level = 100; % �������̖{��
        CCR = unique(CR);
        CCZ = unique(CZ);
        CCR(1) = [];
    %CCZ(1) = [];
    
    
%�ߓ�   
% [RS,ZS,ITYPE,TET,select_phase, Plasma_Current_minus_offset, Psi_z_nega_200, Psi_z_nega_100, Psi_z0, Psi_z100, Psi_z200] = PreUTST(IUTST,WAHAHA,nsemx,eddy_time,SENSio,LWALL,select_phase, time_CCS)
disp(size(reshape(PSI(1:numel(CCR)*numel(CCZ)),numel(CCZ),numel(CCR))))
% figure('Name','Comparison between reconstructed and direct measurement value','NumberTitle','off')
PSI_reconst = 2 * pi * reshape(PSI(1:numel(CCR)*numel(CCZ)),numel(CCZ),numel(CCR));
[Psi_2D, Bz_distribution_plus_EF, Psi_2D_sp]=read_2Darray_data(foldername,shotnum, time_CCS, EF_voltage);

r = [0.11, 0.18, 0.25, 0.32, 0.39, 0.46, 0.53, 0.60, 0.67];
% z = [0.23, 0.18, 0.12, 0.06, 0.00, -0.06, -0.12, -0.18, -0.23];
z = [0.18, 0.12, 0.06, 0.00, -0.06, -0.12, -0.18];
CCZ_idx = round((z + 1) / 0.02 + 1)
CCR_idx = round((r - 0.1) / 0.01 + 1)

fprintf('Error between 2Darray and reconstruction')
Error_2D_reconst =  Psi_2D(2 : 8, 8) - PSI_reconst(CCZ_idx , CCR_idx(8))
Error_2D_reconst_norm = sqrt(sum(Error_2D_reconst.^2))

figure('Name','Comparison between reconstructed and direct measurement value 2','NumberTitle','off')

plot(r(1 : 8), Psi_2D(2, :), 'c--', 'DisplayName', '2D array value at z = 0.18m')
hold on
plot(CCR, PSI_reconst(CCZ_idx(1), :), 'c-', 'DisplayName', 'Reconstructed value at z = 0.18m')
hold on
plot(r(1 : 8), Psi_2D(3, :), 'm--', 'DisplayName', '2D array value at z = 0.12m')
hold on
plot(CCR, PSI_reconst(CCZ_idx(2), :), 'm-', 'DisplayName', 'Reconstructed value at z = 0.12m')
hold on
plot(r(1 : 8), Psi_2D(4, :), 'y--', 'DisplayName', '2D array value at z = 0.06m')
hold on
plot(CCR, PSI_reconst(CCZ_idx(3), :), 'y-', 'DisplayName', 'Reconstructed value at z = 0.06m')
hold on
plot(r(1 : 8), Psi_2D(5, :), 'r--', 'DisplayName', '2D array value at z = 0.00m')
hold on
plot(CCR, PSI_reconst(CCZ_idx(4), :), 'r-', 'DisplayName', 'Reconstructed value at z = 0.00m')
hold on
plot(r(1 : 8), Psi_2D(6, :), 'g--', 'DisplayName', '2D array value at z = -0.06m')
hold on
plot(CCR, PSI_reconst(CCZ_idx(5), :), 'g-', 'DisplayName', 'Reconstructed value at z = -0.06m')
hold on
plot(r(1 : 8), Psi_2D(7, :), 'b--', 'DisplayName', '2D array value at z = -0.12m')
hold on
plot(CCR, PSI_reconst(CCZ_idx(6), :), 'b-', 'DisplayName', 'Reconstructed value at z = -0.12m')
hold on
plot(r(1 : 8), Psi_2D(8, :), 'k--', 'DisplayName', '2D array value at z = -0.18m')
hold on
plot(CCR, PSI_reconst(CCZ_idx(7), :), 'k-', 'DisplayName', 'Reconstructed value at z = -0.18m')
legend

figure('Name','Reconstructed flux on the limiter edge','NumberTitle','off')
plot(CCZ(25 : 76), PSI_reconst(25 : 76 , 2))
hold on
plot(CCZ(25 : 76), PSI_reconst(25 : 76 , 10))
hold on
plot(CCZ(25 : 76), PSI_reconst(25 : 76 , 18))
        
        if 0
           plot(0.50, 0.54, 'ro', 'MarkerFaceColor', 'red')
           hold on
        end

        
        rr = 0.11 : 0.01 : 0.60;
        zz = -0.23 : 0.01 : 0.23
        zz = sort(zz, 'descend')
        
        drawfig = 1
        if drawfig == 1
            figure('Name','Reconstructed flux surface LCFS','NumberTitle','off')    
            dlmwrite(strcat('UTST_',num2str(IUTST),'/Rref.txt'),CCR);
            dlmwrite(strcat('UTST_',num2str(IUTST),'/Zref.txt'),CCZ);
            dlmwrite(strcat('UTST_',num2str(IUTST),'/FLUXref.txt'),reshape(PSI(1:numel(CCR)*numel(CCZ)),numel(CCZ),numel(CCR)));
            PSI_reconst = 2 * pi * reshape(PSI(1:numel(CCR)*numel(CCZ)), numel(CCZ), numel(CCR));
            contourf(CCR , CCZ, PSI_reconst, 300);
            c = colorbar;
            c.Label.String = 'Poloidal flux [mWb]';
            caxis([-0.06 0.06]);
            hold on
            min_PSI_reconst = min(min(PSI_reconst))
            max_PSI_reconst = max(max(PSI_reconst))
            del_contour = (max_PSI_reconst - min_PSI_reconst) / 300
%             contour(CCR, CCZ, PSI_reconst, max(PSI_reconst(35 : 66 , 2)) : del_contour * 0.001 : max(PSI_reconst(35 : 66 , 2)) + del_contour * 0.001, 'r-', 'LineWidth', 2)
            contour(CCR, CCZ, PSI_reconst, max(PSI_reconst(20 : 71 , 2)) : del_contour * 0.001 : max(PSI_reconst(20 : 71 , 2)) + del_contour * 0.001, 'r-', 'LineWidth', 2)
            c = colorbar;
            c.Label.String = 'Poloidal flux [mWb]';
            caxis([-0.06 0.06]);
            hold on
%             contour(rr, zz, Psi_2D_sp, max(Psi_2D_sp(: , 1)) : del_contour * 0.001 : max(Psi_2D_sp(: , 1)) + del_contour * 0.001,...
%             'm-', 'LineWidth', 2)
%             hold on
            contour(rr, zz, Psi_2D_sp, max(PSI_reconst(20 : 71 , 2)) : del_contour * 0.001 : max(PSI_reconst(20 : 71 , 2)) + del_contour * 0.001,...
            'b-', 'LineWidth', 2)   
        end

%%% discrepancy (outboard-side)
%%% search r position PSI_qLCFS on MCCS
points = 0;
points_out = 0;
IdxOut = 1
for i = 1 : 23
    for j = 28 : 61
        if PSI_reconst(i + 39, j) > max(PSI_reconst(20 : 71 , 2))
            continue
        elseif PSI_reconst(i + 39, j) < max(PSI_reconst(20 : 71 , 2))
            delPSI1 = abs(PSI_reconst(i + 39, j - 1) - max(PSI_reconst(20 : 71 , 2)))
            delPSI2 = abs(PSI_reconst(i + 39, j) - max(PSI_reconst(20 : 71 , 2)))
            R_on_reconst(i) = (CCR(j - 1) * delPSI2 + CCR(j) * delPSI1) / (delPSI1 + delPSI2)
            Z_on_reconst(i) = CCZ(i + 39)
            rr_idx = round((CCR(j) - 0.11) / 0.01 + 1)
            
%             error = abs(PSI_reconst(i + 39, j : 51) - Psi_2D_sp(48 - 2 * i, rr_idx : 50))
            % dimensionless (%)
            error = abs((PSI_reconst(i + 39, j : 51) - Psi_2D_sp(48 - 2 * i, rr_idx : 50)) ./ Psi_2D_sp(48 - 2 * i, rr_idx : 50))
            
            Discrepancy(IdxOut) = sum(error)
            points = points + 50 - rr_idx + 1            
            points_out = points_out + 50 - rr_idx + 1            
            IdxOut = IdxOut + 1
            break
        end
    end
end

%%% discrepancy (inboard-side)
%%% search r position PSI_qLCFS on MCCS
points_in = 0
IdxIn = 1
for i = 1 : 23
    for j = 24 : -1 : 2
        if PSI_reconst(i + 39, j) > max(PSI_reconst(20 : 71 , 2))
            continue
        elseif PSI_reconst(i + 39, j) < max(PSI_reconst(20 : 71 , 2))
            delPSI1 = abs(PSI_reconst(i + 39, j + 1) - max(PSI_reconst(20 : 71 , 2)))
            delPSI2 = abs(PSI_reconst(i + 39, j) - max(PSI_reconst(20 : 71 , 2)))
            R_on_reconst(i + 23) = (CCR(j + 1) * delPSI2 + CCR(j) * delPSI1) / (delPSI1 + delPSI2)
            Z_on_reconst(i + 23) = CCZ(i + 39)
            rr_idx = round((CCR(j) - 0.11) / 0.01 + 1)
                        
%             error = abs(PSI_reconst(i + 39, 2 : j) - Psi_2D_sp(48 - 2 * i, 1 : rr_idx))
            % dimensionless (%)
            error = abs((PSI_reconst(i + 39, 2 : j) - Psi_2D_sp(48 - 2 * i, 1 : rr_idx)) ./ Psi_2D_sp(48 - 2 * i, 1 : rr_idx))            
            Discrepancy(IdxOut + 1 + IdxIn) = sum(error)
%             pause
            points = points + rr_idx            
            points_in = points_in + rr_idx
            IdxIn = IdxIn + 1
            break
        end
    end
end


% if drawfig == 1
%     hold on
%     plot(R_on_reconst, Z_on_reconst, 'ro')
%     hold on
% end

for i = 1 : 23
    for j = 35 : 60
        if Psi_2D_sp(2 * i, j) > max(PSI_reconst(20 : 71 , 2))
            continue
        elseif Psi_2D_sp(2 * i, j) < max(PSI_reconst(20 : 71 , 2))
            delPSI1 = abs(Psi_2D_sp(2 * i, j - 1) - max(PSI_reconst(20 : 71 , 2)))
            delPSI2 = abs(Psi_2D_sp(2 * i, j) - max(PSI_reconst(20 : 71 , 2)))
            R_on_2Darray(i) = (rr(j - 1) * delPSI2 + rr(j) * delPSI1) / (delPSI1 + delPSI2)
            Z_on_2Darray(i) = zz(round(2 * i))
            break
        end
    end
end

% if drawfig == 1
%     plot(R_on_2Darray, Z_on_2Darray, 'bo')
%     hold on
% end

        dlmwrite(strcat('UTST_',num2str(IUTST),'/Rrec.txt'),CXK);
        dlmwrite(strcat('UTST_',num2str(IUTST),'/Zrec.txt'),CYK);
        dlmwrite(strcat('UTST_',num2str(IUTST),'/FLUXrec.txt'),SOLK);
        MASK1 = dlmread('MASK1.txt');
        fill(MASK1(:,1),MASK1(:,2),'w')
        hold on
        fill(MASK1(:,1),-MASK1(:,2),'w')
%         xlabel({'r (m)'});
%         ylabel({'z (m)'});
%         axis equal
        if (PLOTW == 0)
            axis([min(VV(:,1)) max(VV(:,1)) -0.5 0.5])%min(VV(:,2)) max(VV(:,2))])
        elseif (PLOTW == 1)
            axis([min(VV(:,1)) max(VV(:,1)) min(VV(:,2)) max(VV(:,2))])%min(VV(:,2)) max(VV(:,2))])
            view(2)
        end
end

axis equal
% discrepancy_out = sum(Discrepancy(1 : 23)) / points_out
% discrepancy_in = sum(Discrepancy(24 : 46)) / points_in
% discrepancy = sum(Discrepancy) / points
discrepancy_out = sum(Discrepancy(1 : IdxOut)) / points_out
discrepancy_in = sum(Discrepancy(IdxOut + 1 : IdxOut + IdxIn)) / points_in
discrepancy = sum(Discrepancy) / points

points
IdxOut
IdxIn
% pause



figure('Name','Determine the position of quasi-LCFS','NumberTitle','off')    
plot(CCR, PSI_reconst(50, :), 'o')
hold on
plot([0, 0.6], [max(PSI_reconst(35 : 66 , 2)), max(PSI_reconst(35 : 66 , 2))], 'b--')
hold on
size_PSI_reconst = size(PSI_reconst);
max(PSI_reconst(35 : 66 , 2))

if 0
    CCR_sub = CCR(25 : 80)
    PSI_reconst_sub = PSI_reconst(50, 25 : 80);
    I = find(PSI_reconst_sub < max(PSI_reconst(35 : 66 , 2)));
    CCR_sub(I(1))
    CCR_sub(I(1) - 1)

    CCZ_sub = CCZ(39 : 67);
    size_PSI_reconst_sub = size(PSI_reconst(39 : 67 , 25 : 80));
    PSI_reconst_sub2 = PSI_reconst(39 : 67 , 25 : 80);
    PSI_q_LCFS = max(PSI_reconst(35 : 66 , 2))
    figure('Name','Determine the position of quasi-LCFS 2','NumberTitle','off')    
    for i = 1 : size_PSI_reconst_sub(1)
        I = find(PSI_reconst_sub2(i, :) < PSI_q_LCFS)
        plot(CCR_sub, PSI_reconst_sub2(i, :), 'ko')
        hold on
        plot(CCR_sub(I(1)), PSI_reconst_sub2(i, I(1)), 'ro')
        hold on
        plot([0.3, 0.8], [max(PSI_reconst(35 : 66 , 2)), max(PSI_reconst(35 : 66 , 2))], 'm-')
        hold on
        if I(1) == 1
            continue
        else
            q_LCFS_reconst_z(i) = CCZ_sub(i);
            r2(i) = CCR_sub(I(1));
            r1(i) = CCR_sub(I(1) - 1);
            PSI_rec_2(i) = PSI_reconst_sub2(i, I(1));
            PSI_rec_1(i) = PSI_reconst_sub2(i, I(1) - 1);
            plot(CCR_sub(I(1) - 1), PSI_reconst_sub2(i, I(1) - 1), 'bo')
            hold on
        end
    end
    rr = rr(20 : 49);
    PSI_q_LCFS_2D = max(Psi_2D_sp(: , 1));
    Psi_2D_sp = Psi_2D_sp(:, 20 : 49);
    size_Psi_2D_sp = size(Psi_2D_sp);
    figure('Name','Determine the position of quasi-LCFS 2 (2Darray)','NumberTitle','off')    
    for i = 1 : size_Psi_2D_sp(1)
        I = find(Psi_2D_sp(i, :) < PSI_q_LCFS_2D)
        plot(rr, Psi_2D_sp(i, :), 'ko')
        hold on
        plot(rr(I(1)), Psi_2D_sp(i, I(1)), 'ro')
        hold on
        plot([0.3, 0.8], [PSI_q_LCFS_2D, PSI_q_LCFS_2D], 'm-')
        hold on
        if I(1) == 1
            continue
        else
            q_LCFS_reconst_2D_z(i) = zz(i)
            r2_2D(i) = rr(I(1));
            r1_2D(i) = rr(I(1) - 1);
            PSI_rec_2D_2(i) = Psi_2D_sp(i, I(1));
            PSI_rec_2D_1(i) = Psi_2D_sp(i, I(1) - 1);
            plot(rr(I(1) - 1), Psi_2D_sp(i, I(1) - 1), 'bo')
            hold on
        end
    %     
    end
    % ((max(PSI_reconst(35 : 66 , 2)) - PSI_rec_1) ./ (PSI_rec_2 - PSI_rec_1))
    q_LCFS_reconst_r = r1 + (r2 - r1) .* ((max(PSI_reconst(35 : 66 , 2)) - PSI_rec_1) ./ (PSI_rec_2 - PSI_rec_1));
    q_LCFS_reconst_2D_r = r1_2D + (r2_2D - r1_2D) .* ((PSI_q_LCFS_2D) - PSI_rec_2D_1) ./ (PSI_rec_2D_2 - PSI_rec_2D_1);
    figure('Name','q-LCFS','NumberTitle','off')    
    plot(q_LCFS_reconst_r, q_LCFS_reconst_z, 'bo')
    hold on
    plot(q_LCFS_reconst_2D_r, q_LCFS_reconst_2D_z, 'ro')

    q_LCFS_reconst_z
    q_LCFS_reconst_r
    q_LCFS_reconst_2D_z

    find_nan = ismissing(q_LCFS_reconst_r)
    QLR_r_idx = find(find_nan == 0)
    QLR_r = q_LCFS_reconst_r(QLR_r_idx)


    q_LCFS_reconst_z(1)
    idx_start = find(q_LCFS_reconst_2D_z > QLR_r(1) - 0.00001 & q_LCFS_reconst_2D_z < QLR_r(1) + 0.00001)

    % idx_start = find(and(q_LCFS_reconst_z < q_LCFS_reconst_2D_z(1) + 0.0001, q_LCFS_reconst_z > q_LCFS_reconst_2D_z(1) - 0.0001))


    % -0.22��Z��0.22m�ŕ]��
    q_LCFS_reconst_z_pick = q_LCFS_reconst_z(1 : 23);
    q_LCFS_reconst_r_pick = q_LCFS_reconst_r(1 : 23);
    q_LCFS_reconst_2D_z_pick = q_LCFS_reconst_2D_z(2 : 2 : 47);
    q_LCFS_reconst_2D_r_pick = q_LCFS_reconst_2D_r(2 : 2 : 47);


    discrepancy = q_LCFS_reconst_r_pick - q_LCFS_reconst_2D_r_pick;


    figure('Name','discrepancy','NumberTitle','off')
    plot(q_LCFS_reconst_z_pick, discrepancy)
    hold on
    plot(q_LCFS_reconst_z_pick, abs(discrepancy))
    hold on
    plot([-0.3, 0.3], [mean(abs(discrepancy)), mean(abs(discrepancy))])

    % calc error area
    size_discrepancy = size(discrepancy)
    for i = 1 : size_discrepancy(2) - 1
        error_area_sub(i) = (abs(discrepancy(i + 1)) + abs(discrepancy(i))) * (q_LCFS_reconst_z_pick(i + 1) - q_LCFS_reconst_z_pick(i)) * 0.5
    end
    error_area = sum(error_area_sub)
end


% q_LCFS_reconst_2D_r = r1_2D + (r2_2D - r1_2D) .* ((PSI_q_LCFS_2D) - PSI_rec_2D_1) ./ (PSI_rec_2D_2 - PSI_rec_2D_1));
% figure('Name','q-LCFS','NumberTitle','off')    
% plot(q_LCFS_reconst_2D_r, q_LCFS_reconst_2D_z)
% 


% C
if (AUTOINPUT(35) == 1)
    AUTO = 0;
end

% if or(AUTO == 0, AUTOINPUT(35) == 0)
%     prompt = '���E��T������? (yes/no) = (1/0)\n';
%     SERX = input(prompt);
% else
%     SERX = AUTOINPUT(36);
% end
% if (SERX == 1)
%     fid50 = fopen(strcat('UTST_',num2str(IUTST),'/BOUNDARY_P1.TXT'),'w'); %50
%     fid51 = fopen(strcat('UTST_',num2str(IUTST),'/BOUNDARY_P2.TXT'),'w'); % 51
% ZERO=0.0;
% for J = 1:1000
%     fprintf(fid51,myformat,ZERO,ZERO);%  ! ##1
%     fprintf(fid50,myformat,ZERO,ZERO);% ! ##1
% end
% [statR,endR,AAA1,BBB1,PPPSI1,LOCAT,CCSNEWTON] = SADDLE(CR,CZ,PSI,Z0,RR,FFOUT,...
% n100,n50,RCCS,ZCCS,RC,ZC,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,...
% ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,GETA,WAHAHA,MXCCS,RMYU0,NE,NCCS,fid99,fid100,CCS,AUTO,AUTOINPUT);
% %
% %
% [TOTMOD,TRUTOT,EPSAX,DEFR,DEFZ,C0R,C0Z,I1I0,CAPPER] = NEWTON(IUTST,ICONT,RR,RCCS,ZCCS,...
% FFOUT,n100,n50,R0,Z0,RC,ZC,ECIGRP,ECI,KCMX,AAA1,BBB1,...
% PPPSI1,TOTMOD,I1I0,RCCN,ZCCN,NCCN,SOU,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,...
% AMYU0,NONC,GETA,LOCAT,NCCS,RMYU0,fid99,fid100,fid50,fid51,WAHAHA,MXCCS,CCS,NE,CCSNEWTON,AUTO,AUTOINPUT);
% %
% else
%     CCSNEWTON  = 0;
% end
% 
% 
% 
% % PSI_reconst = reshape(PSI(1:numel(CCR)*numel(CCZ)),numel(CCZ),numel(CCR));
% % plot(PSI_reconst(70, :))
% % 
% 
% 
% 
% prompt ='Plot the LCMS or plasma boundary? (Yes/No)=(1/0)\n';
% PLOTB = input(prompt);
%     if (PLOTB > 0)
%         if and(IUTST >= 4, IUTST ~= 10)
%             BOUND_ans = dlmread(strcat('UTST_',num2str(IUTST),'/PS.txt'));
%             hold on
% %             plot(BOUND_ans(:,1),BOUND_ans(:,2),'-b','LineWidth',2);
%         elseif (IUTST == 10)
%             BOUND_ans = dlmread(strcat('UTST_',num2str(IUTST),'/PS.txt'));
% %             hold on
% %             plot(vertcat(BOUND_ans(1:216,1),BOUND_ans(1,1)),vertcat(BOUND_ans(1:216,2),BOUND_ans(1,2)),'-b','LineWidth',2);
% %             hold on
% %             plot(vertcat(BOUND_ans(217:end,1),BOUND_ans(217,1)),vertcat(BOUND_ans(217:end,2),BOUND_ans(217,2)),'-b','LineWidth',2);
%         end
%         BOUND = dlmread(strcat('UTST_',num2str(IUTST),'/BOUNDARY_P',num2str(CCSNEWTON),'.TXT'));
%         hold on
%         plot(BOUND(:,1),BOUND(:,2),'-r','LineWidth',2);
%         prompt = 'Make the data for plasma current density reconstruction? (Yes/No)=(1/0)\n';
%         DATA_P = input(prompt);
%         if (DATA_P > 0)
%             MakePlasmaCurrentDensityData(GETA,RCCS,ZCCS,FFOUT,NCCN,NE,CCS,BOUND,MXCCS,IUTST,n100,...
%             n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,...
%             fid99,fid100,RMYU0)
%         end
%     end
% prompt ='Again the LCMS Search (0)?  or Quit (1)?\n';
% LTRY = input(prompt);
% while (LTRY == 0)
%     [statR,endR,AAA1,BBB1,PPPSI1,LOCAT,CCSNEWTON] = SADDLE(CR,CZ,PSI,Z0,RR,FFOUT,...
%     n100,n50,RCCS,ZCCS,RC,ZC,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,...
%     ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,GETA,WAHAHA,MXCCS,RMYU0,NE,NCCS,fid99,fid100,CCS,AUTO,AUTOINPUT);
% % 
% %
%     [TOTMOD,TRUTOT,EPSAX,DEFR,DEFZ,C0R,C0Z,I1I0,CAPPER] = NEWTON(IUTST,ICONT,RR,RCCS,ZCCS,...
%     FFOUT,n100,n50,R0,Z0,RC,ZC,ECIGRP,ECI,KCMX,AAA1,BBB1,...
%     PPPSI1,TOTMOD,I1I0,RCCN,ZCCN,NCCN,SOU,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,...
%     AMYU0,NONC,GETA,LOCAT,NCCS,RMYU0,fid99,fid100,fid50,fid51,WAHAHA,MXCCS,CCS,NE,CCSNEWTON,AUTO,AUTOINPUT);
% % 
% %  
% prompt ='Plot the LCMS or plasma boundary? (Yes/No)=(1/0)\n';
% PLOTB = input(prompt);
%     if (PLOTB > 0)
%         if (CCSNEWTON > 0)
%             if (IUTST >= 4)
%                 BOUND_ans = dlmread(strcat('UTST_',num2str(IUTST),'/PS.txt'));
% %                 hold on
% %                 plot(BOUND_ans(:,1),BOUND_ans(:,2),'-b','LineWidth',2);
%             elseif (IUTST == 10)
%                 BOUND_ans = dlmread(strcat('UTST_',num2str(IUTST),'/PS.txt'));
%                 hold on
%                 plot(vertcat(BOUND_ans(1:216,1),BOUND_ans(1,1)),vertcat(BOUND_ans(1:216,2),BOUND_ans(1,2)),'-b','LineWidth',2);
%                 hold on
%                 plot(vertcat(BOUND_ans(217:end,1),BOUND_ans(217,1)),vertcat(BOUND_ans(217:end,2),BOUND_ans(217,2)),'-b','LineWidth',2);
%             end
%             BOUND = dlmread(strcat('UTST_',num2str(IUTST),'/BOUNDARY_P',num2str(CCSNEWTON),'.TXT'));
%             hold on
%             plot(BOUND(:,1),BOUND(:,2),'-r','LineWidth',2);
%         else
%             BOUND_ans1 = dlmread(strcat('ST40_',num2str(IST40),'/@ST40_boundary1.txt'));
%             hold on
%             plot(BOUND_ans1(:,1),BOUND_ans1(:,2),'-b','LineWidth',2);
%             BOUND1 = dlmread(strcat('UTST_',num2str(IUTST),'/BOUNDARY_P1.TXT'));
%             hold on
%             plot(BOUND1(:,1),BOUND1(:,2),'-r','LineWidth',2);
%             BOUND2 = dlmread(strcat('UTST_',num2str(IUTST),'/BOUNDARY_P2.TXT'));
%             hold on
%             plot(BOUND2(:,1),BOUND2(:,2),'-r','LineWidth',2);
%         end
%     end
%     prompt ='Again the LCMS Search (0)?  or Quit (1)?\n';
%     LTRY = input(prompt);
% end
% prompt ='Plot the error around boundary? (Yes/No)=(1/0)\n';
% PERROR = input(prompt);
% if (PERROR > 0)
%     ERROR_BOUNDARY(IUTST,GETA,RCCS,ZCCS,FFOUT,n100,n50,RC,ZC,...
%     ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,...
%     RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS,MXCCS);
% end
% % C      
% prompt ='Modify the CCS shape (0)?  or Quit (1)?\n';
% MODCCS = input(prompt);
% if (MODCCS == 0)
%     R0 = C0R;
%     Z0 = C0Z;
%     ICONT = ICONT+1;
%     GHR(1:NCCS) = DEFR(1:NCCS);
%     GHZ(1:NCCS) = DEFZ(1:NCCS);
%     EPS0 = 0.001;
%     fprintf('%s%d\r\n','EPSAX=',EPSAX)
% 	if (EPSAX < EPS0)
%         fprintf('%s\r\n','Congraturation!!!!!');
%         fprintf('%s\r\n','CCS Modification Search is Over!!');
% %        8888 CALL CPU_TIME(CTM)
% %        WRITE(*,'(''***CPU Time ='',1PE11.4)') CTM
%  	else
% % 	      GOTO 9999  ! CCS����蒼����FORM����ĊJ�B
%  	end
% else
% end
% fclose('all');
end
%
%
%% D_CCS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111 
function[RCCS,ZCCS,RCCN,ZCCN,NCCN,NE,NCCS] = D_CCS(R0,Z0,RR,CAPPER,MXCCS,CCS,MSEC)
RCCS = zeros(CCS,MXCCS+1);
ZCCS = zeros(CCS,MXCCS+1);%1
RCCN = zeros(CCS,MXCCS+1);
ZCCN = zeros(CCS,MXCCS+1);
SR(1:CCS,1) = R0-RR;
SR(1:CCS,3) = SR(1:CCS,1);
SR(1:CCS,2) = R0+RR;
SZ(1:CCS,1)=Z0+RR.*CAPPER;
SZ(1:CCS,2)=Z0;
SZ(1:CCS,3)=Z0-RR.*CAPPER;
NE = MSEC(1:CCS,1)+MSEC(1:CCS,2)+MSEC(1:CCS,3);
NCCS = NE*2;
fid58 = fopen('@D_CCS_SRpoints.txt','w');
fid59 = fopen('@D_CCS_CheckWrite.txt','w');
fid13 = fopen('DiscontinuousNodePoints.txt','w');
frewind(fid13);
myformat = '%d %d\r\n';
for III = 1:CCS
    fprintf(fid58,myformat,R0(III),Z0(III));
    for I = 1:3
        fprintf(fid58,myformat, SR(III,I),SZ(III,I));
    end
    %MSEC(1)=1; MSEC(2)=1; MSEC(3)=1; % �R�ɕ�����ꂽ�Z�N�V�����̂��ꂼ��̋��E�v�f��
    %MSEC(1)=2; MSEC(2)=2; MSEC(3)=2;
    %
    %% �Ȑ����̃��b�V���_���`
    II=0;
    for I = 1:2;
        M2 = MSEC(III,I)*2; % ���b�V���_�̐�
        DEL = 1.0/M2;
        GISTAT = -1.0+(I-1); % �O�U�C -1 0 1
        for J = 1:M2
            GI = GISTAT+DEL*(J-1); %�O�U�C �����E�v�f���ɕ������Ă���
            % Compute the values of the shape functions at the integration points
            F1 = GI*(GI-1.0)/2;
            F2 = 1.0-GI^2;
            F3 = GI*(GI+1.0)/2;
            II = II+1;
            RCCS(III,II) = SR(III,1)*F1+SR(III,2)*F2+SR(III,3)*F3;
            ZCCS(III,II) = SZ(III,1)*F1+SZ(III,2)*F2+SZ(III,3)*F3; 
        end          
    end 
    %% �������̃��b�V���_���`
    II=II+1;
    RCCS(III,II)=SR(III,3);
    ZCCS(III,II)=SZ(III,3); 
    KADO=II;
    M2=MSEC(III,3)*2;
    DELR=(SR(III,1)-SR(III,3))/M2;
    DELZ=(SZ(III,1)-SZ(III,3))/M2;
    for J = 1:M2
        II=II+1;
        RCCS(III,II)=RCCS(III,II-1)+DELR;
        ZCCS(III,II)=ZCCS(III,II-1)+DELZ;
    end
    NCCS(III) = II-1;
    NCCN(III) = NE(III)*3;
    fprintf('%d%s%d  %d  %d\r\n',III,'�Ԗ�CCS: NE/NCCS/NCCN = ',NE(III), NCCS(III),NCCN(III));
    for II=1:NCCS(III)+1
        fprintf('%d %d %d\r\n',II,RCCS(III,II),ZCCS(III,II));
        fprintf(fid59, myformat, RCCS(III,II), ZCCS(III,II));
        if (II == KADO)
 %           fprintf(fid59,'%s\r\n', '==');
            fprintf(fid59,myformat,RCCS(III,II), ZCCS(III,II));
        else
        end
    end
    %% CCS��̔�K���v�f�ߓ_���W�̍쐬
    RCCS(III,NCCS+1) = RCCS(III,1);
    ZCCS(III,NCCS+1) = ZCCS(III,1);
    NCCN(III) = NE(III)*3; % ��K���v�f
    %
    I = 1:NE(III);%NE=��K���v�f��
    RCCN(III,3*I-2)=(5*RCCS(III,2*I-1)+5.*RCCS(III,2*I)-RCCS(III,2*I+1))/9;
    ZCCN(III,3*I-2)=(5*ZCCS(III,2*I-1)+5.*ZCCS(III,2*I)-ZCCS(III,2*I+1))/9;
    RCCN(III,3*I-1)=RCCS(III,2*I);
    ZCCN(III,3*I-1)=ZCCS(III,2*I);
    RCCN(III,3*I)  =(5*RCCS(III,2*I+1)+5.*RCCS(III,2*I)-RCCS(III,2*I-1))/9;
    ZCCN(III,3*I)  =(5*ZCCS(III,2*I+1)+5.*ZCCS(III,2*I)-ZCCS(III,2*I-1))/9;
%
    for I = 1:NCCN(III);
        fprintf(fid13,myformat,RCCN(III,I), ZCCN(III,I));
    end;
end
fclose(fid13);
fclose(fid58);
fclose(fid59);
end
%
%
% Psi_z��GM���o�͂��������̌a�������z�Ƀt���b�N�X���[�v�ɂ����������l�ŕ␳�������́i�����l��^�����j
function [RS,ZS,ITYPE,TET,select_phase, Plasma_Current_minus_offset, Psi_z_nega_200, Psi_z_nega_100, Psi_z0, Psi_z100, Psi_z200, R_2d] = PreUTST(IUTST,WAHAHA,nsemx,eddy_time,SENSio,LWALL,select_phase, time_CCS, foldername_TF, shotnum_TF, foldername, shotnum, EF_voltage)     
%function [RS,ZS,ITYPE,TET] = PreUTST(IUTST,WAHAHA,nsemx,LWALL,IPPP,L0BR)     
IMAX = 512;   % IMAX=513 (�d���f�[�^�폜)
MAXM = 10;
RSENS = 0.113150; 
RWALL = 0.108150;  
% ***********************************************************************
% �s��̑傫�����w�肵�Ȃ����������\��������̂ŃR�����g�A�E�g
%RS = zeros(1,nsemx); 
%ZS = zeros(1,nsemx); 
%ITYPE = zeros(1,nsemx); 
%TET = zeros(1,nsemx);
% ***********************************************************************
%     =================================================================
% probe & flux-loop data(position data)
%     -----------------------------------------------------------------
% NTPB   : number of T-probes.
% NNPB   : number of N-probes.
% NAPB   : number of all(T & N) probes.
% NFLX   : number of Flux-Loops.
%     =================================================================
TPRB = zeros(1,1024); 
FLXLP = zeros(1,1024);
%
% UTST �̕ǂ̊􉽍\��������
RSEC = [0.694D0, 0.694D0, 0.5985D0, 0.5985D0, 0.10815D0, 0.10815D0,... 
        0.10815D0, 0.5985D0, 0.5985D0, 0.694D0,  0.694D0];
ZSEC = [0.0D0, 0.285D0, 0.285D0, 0.9985D0, 0.9985D0, 0.0D0,...
       -0.9985D0, -0.9985D0, -0.285D0, -0.285D0 , 0.0D0];    
   
   
% foldername_TF = '180528';
% shotnum_TF = '001';
% foldername = '180913';
% shotnum = '001';
[time, Plasma_Current_TFshot, Coil_Current_TFshot, Flux_Loop_TFshot, Magnetic_Probe_TFshot, Magnetic_Probe_Lowpass_TFshot] = read_CCS_data(foldername_TF, shotnum_TF);   
% foldername = '180515';
% shotnum = '009';
% foldername = '180913';
% shotnum = '007';
[time, Plasma_Current, Coil_Current, Flux_Loop, Magnetic_Probe, Magnetic_Probe_Lowpass] = read_CCS_data(foldername, shotnum);   


size(Plasma_Current)

Plasma_Current_minus_offset = Plasma_Current - mean(Plasma_Current(1000 : 2000));
% plot(time, Plasma_Current)
% hold on
% 100us�ňړ�����
% figure('Name','Plasma Current smooth','NumberTitle','off')
Plasma_Current_minus_offset_smooth = smoothdata(Plasma_Current_minus_offset, 'gaussian', 200);
% plot(time, Plasma_Current_minus_offset_smooth)
% xlim([8.5, 10.5])
% 

% figure('Name','Plasma Current','NumberTitle','off')
% subplot(5, 1, 5); plot(time, Plasma_Current_minus_offset, 'color', [0.7 0.7 0.7], 'DisplayName', 'Plasma current (raw)')
% hold on
% subplot(5, 1, 5);plot(time, Plasma_Current_minus_offset_smooth, 'r', 'DisplayName', 'Plasma current (smoothing)')
% xlim([8.5, 10.5])
% ylim([-100, 150])
% xlabel('Time [ms]');
% ylabel('Current [kA]');
% legend({}, 'FontSize', 5, 'Location', 'eastoutside')

% hold on
% plot(time_CCS, Plasma_Current_minus_offset(round(time_CCS / (0.5 * 0.001))), 'o')

PF_Current_for_CCS = zeros(4, 1);
% figure('Name','Coil Current','NumberTitle','off')
% smoothdata(Plasma_Current_minus_offset, 'gaussian', 200, 'DisplayName', 'Plasma Current smoothed')
for i = 1 : 4
    % �㉺�̃R�C���̓d���l�̕��ς�p����
    PF_Current_for_CCS(i, 1) = (Coil_Current(2 * i - 1, round(time_CCS / (0.5 * 0.001))) + Coil_Current(2 * i, round(time_CCS / (0.5 * 0.001)))) / 2;
    % �ۂ�
%     PF_Current_for_CCS(i, 1) = round(PF_Current_for_CCS(i, 1), 5);
%     subplot(4, 1, i); plot(time, Coil_Current(2 * i - 1, :), time_CCS, Coil_Current(2 * i - 1, round(time_CCS / (0.5 * 0.001))), 'o', 'DisplayName', 'PF%d (upper side)')
    subplot(5, 1, i); plot(time, Coil_Current(2 * i - 1, :), 'color', [0.7 0.7 0.7], 'DisplayName', 'PF#   current (upper side)')
    xlim([8.5, 10.5])
    ylim([-50, 100])
    hold on
    subplot(5, 1, i); plot(time, Coil_Current(2 * i, :), 'color', [0.9 0.9 0.9], 'DisplayName', 'PF#   current (lower side)')
    hold on
    subplot(5, 1, i); plot(time, smoothdata(Coil_Current(2 * i - 1, :), 'gaussian', 200), 'r', 'DisplayName', 'Smoothed PF#   current (upper side)')
    hold on
    subplot(5, 1, i); plot(time, smoothdata(Coil_Current(2 * i, :), 'gaussian', 200), 'b', 'DisplayName', 'Smoothed PF#   current (lower side)')
    hold on
    legend({}, 'FontSize', 5, 'Location', 'eastoutside')
    
end
fprintf('%d', PF_Current_for_CCS)



MP_size = size(Magnetic_Probe);
BZ_t = zeros(MP_size(1), MP_size(2));
BZ_t_TFshot = zeros(MP_size(1), MP_size(2));
BZ_t_Lowpass = zeros(MP_size(1), MP_size(2));
BZ_t_minusTF = zeros(MP_size(1), MP_size(2));



fileid = fopen('TF.txt', 'r');
TF = fscanf(fileid, '%f');

fileid2 = fopen('Bz_EF_at_MP_clockwise.txt', 'r');
EF_MP = fscanf(fileid2, '%f');

fileid3 = fopen('psi_EF_at_FL_clockwise.txt', 'r');
EF_FL = fscanf(fileid3, '%f');


% �Z���T�ʒu�ł�EF�̊�^���v�Z
fid61 = fopen('CCS_FLXLP_sensor_position.txt','r');
fid68 = fopen('CCS_MP_sensor_position.txt','r');
textscan(fid61,'%s',1,'Delimiter','\n');
z_flxlp = zeros(1, 35);
r_flxlp = ones(1, 35);
for i = 1:35
    temp = textscan(fid61,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    r_flxlp(1, i) = str2double(temp_m(1));
    z_flxlp(1, i) = str2double(temp_m(2));
end
% �����t���b�N�X���[�v��r=0.108785m�ɐݒu
r_flxlp = r_flxlp .* 0.108785;

z_mp = zeros(1, 40);
r_mp = ones(1, 40);
textscan(fid68,'%s',1,'Delimiter','\n');
for i = 1:40
    temp = textscan(fid68,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    r_mp(1, i) = str2double(temp_m(1));
    z_mp(1, i) = str2double(temp_m(2));
end

[Bz_EF_at_sensor_f, Psi_EF_at_sensor_f] = EF_calc_for_CCS_probe(r_flxlp, z_flxlp, 1, 19, EF_voltage);
[Bz_EF_at_sensor_b, Psi_EF_at_sensor_b] = EF_calc_for_CCS_probe(r_mp, z_mp, 1, 40, EF_voltage);


% plot(EF_MP)
% 
% plot(EF_FL)
% 
% plot(TF(1 : 30000))
% 

TF_txt = zeros(40, 30000);
BZ_Distribution = zeros(1, MP_size(1));
EF_Distribution_MP = zeros(1, MP_size(1));
figure('Name','Axial Magnetic Field','NumberTitle','off')
for i = 1 : MP_size(1)
    TF_txt(i, :) = TF(1 + (i - 1) * 30000 : i * 30000) * 0.001;
    BZ_t(i, :) = cumtrapz(time, Magnetic_Probe(i, :));
    BZ_t_Lowpass(i, :) = cumtrapz(Magnetic_Probe_Lowpass(i, :) - mean(Magnetic_Probe_Lowpass(i, 1:1000))) * 0.5e-6 * 242;
%     BZ_t_TFshot(i, :) = cumtrapz(Magnetic_Probe_TFshot(i, :) - mean(Magnetic_Probe_TFshot(i, 1:1000))) * 0.5e-6 * 242;
    BZ_t_TFshot(i, :) = cumtrapz(Magnetic_Probe_Lowpass_TFshot(i, :) - mean(Magnetic_Probe_Lowpass_TFshot(i, 1:1000))) * 0.5e-6 * 242;
%     BZ_t_TFshot(i, :) = cumtrapz(Magnetic_Probe_TFshot(i, :)) * 0.5e-6 * 242;
%     BZ_t_minusTF(i, :) = BZ_t_Lowpass(i, :) - BZ_t_TFshot(i, :);
    BZ_t_minusTF(i, :) = BZ_t_Lowpass(i, :) - BZ_t_TFshot(i, :);    
%     subplot(5, 8, i ); plot(time, BZ_t(i, :), 'DisplayName', 'BZ + BT')
%     hold on
    subplot(5, 8, i ); plot(time, BZ_t_Lowpass(i, :) - BZ_t_TFshot(i, :), 'LineWidth', 2, 'DisplayName', 'BZ + BT (lowpass) - TF')
    hold on
    subplot(5, 8, i ); plot(time, BZ_t_TFshot(i, :), 'DisplayName', 'BT')
    hold on
    subplot(5, 8, i ); plot(time, TF(1 + (i - 1) * 30000 : i * 30000) * 0.001, 'DisplayName', 'BT txt')
    hold on
%     subplot(5, 8, i ); plot(time, BZ_t(i, :) - BZ_t_TFshot(i, :), 'DisplayName', 'BZ')
%     hold on
    subplot(5, 8, i ); plot(time, BZ_t_minusTF(i, :), 'DisplayName', 'BZ (lowpass) - BT_txt')
    subplot(5, 8, i ); plot(time_CCS, BZ_t_minusTF(i, round(time_CCS / (0.5 * 0.001))), 'o')
    BZ_Distribution(1, i) = BZ_t_minusTF(i, round(time_CCS / (0.5 * 0.001)));
%     EF_Distribution_MP(1, i) = EF_MP(i);
%     Bz_EF_at_sensor_b
%     EF_Distribution_MP(1, i) = Bz_EF_at_sensor_b(1, i);    
    if i == 1
        xlabel('Time [us]');
        ylabel('Axial Magnetic Field [T]');
        legend
    end
end

size(BZ_Distribution)
size(EF_MP)
size(Bz_EF_at_sensor_b)
EF_MP
Bz_EF_at_sensor_b

BZ_Distribution = BZ_Distribution + Bz_EF_at_sensor_b;
% BZ_Distribution = BZ_Distribution + EF_Distribution_MP;

figure('Name','Axial Field Distribution (all) including EF effect','NumberTitle','off')
BZ_Distribution(30) = [];
plot(BZ_Distribution)


FL_size = size(Flux_Loop);
PSI_t = zeros(FL_size(1), FL_size(2));
PSI_Distribution = zeros(1, FL_size(1));
EF_Distribution_FL = zeros(1, FL_size(1));
figure('Name','Poloidal Flux','NumberTitle','off')
% FL signal [V] = Gain * 1 / (RC) * (integral((d��/dt)dt))
% Gain : 100 (inboard), 10 (outboard)
% RC : 1 (inboard), 470 (outboard)
for i = 1 : FL_size(1)
    if and(1 <= i, i <= 19)
%         PSI_t(i, :) = (510. * 0.001 / 100.) * Flux_Loop(i, :);
        
        % �t���b�N�X���[�v��TF��������190210�i�ߓ��j1
        PSI_t(i, :) = (Flux_Loop(i, :)  - Flux_Loop_TFshot(i, :)) / 202;
%         PSI_t(i, :) = Flux_Loop(i, :) / 202;
        PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
        subplot(5, 7, i ); plot(time, PSI_t(i, :), 'DisplayName', 'inboard')
        legend
        PSI_Distribution(1, i) = PSI_t(i, round(time_CCS / (0.5 * 0.001)));
       EF_Distribution_FL(1, i) = EF_FL(i);    
    else
%         PSI_t(i, :) = (470. * 0.001 / 10.) * Flux_Loop(i, :);
%         PSI_t(i, :) = Flux_Loop(i, :) / 37;
        % �t���b�N�X���[�v��TF��������190210�i�ߓ��j
        PSI_t(i, :) = (Flux_Loop(i, :)  - Flux_Loop_TFshot(i, :)) / 37;
        PSI_t(i, :) = PSI_t(i, :) - mean(PSI_t(i, 1000:2000));
        subplot(5, 7, i ); plot(time, PSI_t(i, :), 'DisplayName', 'outboard, top and bottom')
        PSI_Distribution(1, i) = PSI_t(i, round(time_CCS / (0.5 * 0.001)));
%         EF_Distribution_FL(1, i) = EF_FL(i);    
%         EF_Distribution_FL(1, i) = Psi_EF_at_sensor_f(i);    
    end
end

% PSI_Distribution = PSI_Distribution + EF_Distribution_FL;
PSI_Distribution = PSI_Distribution + Psi_EF_at_sensor_f;


figure('Name','Poloidal Flux Distribution (all) including EF effect','NumberTitle','off')
PSI_Distribution(6) = [];
PSI_Distribution(35 - 1) = [];
plot(PSI_Distribution)


   
   
   
   
   

if (select_phase == 0)
%     fid61 = fopen('Parameters_FL_crockwise_180515010_t9500_inpsizero.txt','r');
%     fid68 = fopen('Parameters_MP_crockwise_180515010_t9500_inpsizero.txt','r');
    fid61 = fopen('Parameters_FL_crockwise_180515010_t9500.txt','r');
    fid68 = fopen('Parameters_MP_crockwise_180515010_t9500.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
%    IMAX_f = 1;
elseif (select_phase == 1)
%     fid61 = fopen('Parameters_FL_clockwise_180515010_t9650.txt','r');
%     fid68 = fopen('Parameters_MP_clockwise_180515010_t9650.txt','r');

%inboard�����l���X���[�W���O��������
%     fid61 = fopen('Parameters_FL_clockwise_180515010_t9650_inboard_smoothing.txt','r');
%     fid68 = fopen('Parameters_MP_clockwise_180515010_t9650.txt','r');
%     
%     fid61 = fopen('Parameters_FL_clockwise_180515010_t9650_2.txt','r');
%     fid68 = fopen('Parameters_MP_clockwise_180515010_t9650_2.txt','r');
%     fid61 = fopen('Parameters_FL_crockwise_180515010_t9650_inpsizero.txt','r');
%     fid68 = fopen('Parameters_MP_crockwise_180515010_t9650_inpsizero.txt','r');
    fid61 = fopen('Parameters_FL_clockwise_180515010_t9650_EFkaizen.txt','r');
    fid68 = fopen('Parameters_MP_clockwise_180515010_t9650_EFkaizen.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
%    IMAX_f = 10;
elseif (select_phase == 2)
    fid61 = fopen('Parameters_FL_crockwise_180515010_t9550.txt','r');
    fid68 = fopen('Parameters_MP_crockwise_180515010_t9550.txt','r');
%     fid61 = fopen('Parameters_FL_crockwise_180515010_t9650_inpsizero.txt','r');
%     fid68 = fopen('Parameters_MP_crockwise_180515010_t9650_inpsizero.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
%    IMAX_f = 10;
elseif (select_phase == 3)
    fid61 = fopen('Parameters_FL_clockwise_180515010_t9450_EFkaizen.txt','r');
    fid68 = fopen('Parameters_MP_clockwise_180515010_t9450_EFkaizen.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
%    IMAX_f = 10;
elseif (select_phase == 4)
    fid61 = fopen('Parameters_FL_clockwise_180515010_t9530_EFkaizen.txt','r');
    fid68 = fopen('Parameters_MP_clockwise_180515010_t9530_EFkaizen.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
elseif (select_phase == 5)
%     fid61 = fopen('Parameters_FL_clockwise_180515010_t8000.txt','r');
%     fid68 = fopen('Parameters_MP_clockwise_180515010_t8000.txt','r');
    fid61 = fopen('Parameters_FL_clockwise_180515010_t8000_EFkaizen.txt','r');
    fid68 = fopen('Parameters_MP_clockwise_180515010_t8000_EFkaizen.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
elseif (select_phase == 6)
    fid61 = fopen('Parameters_FL_clockwise_180515010_t8000_EFkaizen.txt','r');
    fid68 = fopen('Parameters_MP_clockwise_180515010_t8000_EFkaizen.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
elseif (select_phase == 7)
    fid61 = fopen('Parameters_FL_clockwise_180515010_t8000_EFkaizen.txt','r');
    fid68 = fopen('Parameters_MP_clockwise_180515010_t8000_EFkaizen.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
elseif (select_phase == 8)
    fid61 = fopen('Parameters_FL_clockwise_180515010_t8000_EFkaizen.txt','r');
    fid68 = fopen('Parameters_MP_clockwise_180515010_t8000_EFkaizen.txt','r');
    IMAX_b = 39;
    IMAX_f = 33;
end




fid60 = fopen('@UTST_CheckWrite.txt','w');
%%%����fid62���v���n�̍��W�ɂ���
fid62 = fopen('@UTST_SenPos.txt','w');
fid63 = fopen('@UTST_VECTOR.txt','w');
fid64 = fopen('@UTST_WallGeom.txt','w');
fid65 = fopen('@UTST_CoilGeom.txt','w');
% ���C�f�[�^�ǂݍ���
textscan(fid61,'%s',1,'Delimiter','\n');  % ��s�Ƃ΂�
%for I = 1:IMAX
for I = 1:IMAX_f
    temp = textscan(fid61,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});   
    R = str2double(temp_m(1));
    Z = str2double(temp_m(2));
    fprintf(fid62,'%d %d\r\n',R,Z);
    fprintf(fid62,'%d %d\r\n',R,-Z);
end
%
for I = 1:MAXM+1
    fprintf(fid64,'%d %d\r\n',RSEC(I),ZSEC(I));
end
%
fprintf(fid65,'%d %d\r\n', 0.80, +1.07); 
fprintf(fid65,'%d %d\r\n', 0.80, -1.07);      
% ***************************************************** %
%   Flux Loop                                           %
% ***************************************************** %
frewind(fid61);
textscan(fid61,'%s',1,'Delimiter','\n');  %�������傤�Ƃ΂�
II = 0;
%
LGP = 1; % �����f�[�^���΂��������Ō���I�I�I�I17 15
LGPH = floor(LGP/2)+1;
LG = 0;

for I = 1:IMAX_f
    temp = textscan(fid61,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    R = str2double(temp_m(1));
    Z = str2double(temp_m(2));
%     PSI = str2double(temp_m(3));    
    PSI = PSI_Distribution(1, I);
    
    LG = LG+1;
    LLG = LG;
    if (LG == LGP) % 1�����Ń��Z�b�g����
        LG = 0;
    end
    if (LLG ~= LGPH) % 1������܂ł͓ǂ݂Ƃ΂�
        continue  
    end
     if and(LWALL > 0, and (R < 0.12, abs(Z) < 0.9)); % �ǉ����̃Z���T�[����2016.9.5ushiki
         continue
     end
    % �s�K�v�ȃt���b�N�X���[�v��p��
    
%     if R < 0.2
%         continue
%     end
    
    % inbord���̃t���b�N�X���[�v��p��
%     if and(R < 0.13, abs(Z) < 1) % �ǉ����̃Z���T�[����
%          continue
%     end    
%     if or(and(R == 0.689126, abs(Z) < 0.085),and(and(R ~= 0.689126, R ~= 0.113024), abs(Z) < 0.25)); % �ǉ����̃Z���T�[����
%          continue
%     end
    E = abs((R-RSENS)/RSENS);
%     if and(LWALL > 0, E < 1.0D-5)
%         II = II+1;
%         % (IPPP == 0) �ǂ�PSI��0�ɂ��� (IPPP ~= 0) �ʏ�
%             RS(2*II-1)=RWALL.*(IPPP == 0) + R.*(IPPP ~= 0);
%             RS(2*II)=RWALL.*(IPPP == 0) + R.*(IPPP ~= 0);
%             FLXLP(2*II-1)=0.0D0.*(IPPP == 0) + PSI.*(IPPP ~= 0);
%             FLXLP(2*II)=0.0D0.*(IPPP == 0) + PSI.*(IPPP ~= 0);
%     else
        II=II+1;
        if (IUTST <= 9)
            RS(2*II-1)=R;
            RS(2*II)=R;
            FLXLP(2*II-1)=PSI;
            FLXLP(2*II)=PSI;
%     end
            ZS(2*II-1)=Z;
            ZS(2*II)=-Z;
            TET(2*II-1)=0.0D0;
	        TET(2*II)=0.0D0;
	        ITYPE(2*II-1)=0;
	        ITYPE(2*II)=0;
	        fprintf(fid60,'%d %d %d %d\r\n', 2*II-1,RS(2*II-1),ZS(2*II-1),FLXLP(2*II-1));
            fprintf(fid60,'%d %d %d %d\r\n',2*II,RS(2*II),ZS(2*II),FLXLP(2*II));
        elseif (IUTST == 10)
            RS(II) = R;
            FLXLP(II) = PSI;
%             FLXLP(II) = PSI / (2 * pi);
            ZS(II) = Z;
            TET(II) = 0.0D0;
	        ITYPE(II) = 0;
	        fprintf(fid60,'%d %d %d %d\r\n', II,RS(II),ZS(II),FLXLP(II));
        elseif (IUTST == 11)
            RS(II) = R;
            FLXLP(II) = PSI;
            ZS(II) = Z;
            TET(II) = 0.0D0;
	        ITYPE(II) = 0;
	        fprintf(fid60,'%d %d %d %d\r\n', II,RS(II),ZS(II),FLXLP(II));
        end
end% 100

% �ߓ�
z_sp_max = 0.6;
z_sp_min = -0.6;
del_z_sp = 0.01;
zz = z_sp_min : del_z_sp : z_sp_max; % zz = 0.01 * (n - 1) + (-0.60)
% 2�����A���C�f�[�^�ɑ������킹�鎥���𒊏o�@���@FLXLP_inboard_spline
% zz_pickup = 0.0;
zz_pickup = [-0.2, -0.1, 0, 0.1, 0.2];
zz_idx = round((zz_pickup + 0.60) / 0.01 + 1);
FLXLP_inboard_spline = spline(ZS(1 : 18), FLXLP(1 : 18), zz);
figure('Name','Axial Distribution Inboard-side Flux','NumberTitle','off')
plot(ZS(1 : 18), FLXLP(1 : 18), 'ko')
hold on
plot(zz, FLXLP_inboard_spline, 'r-')
hold on
plot(zz(zz_idx), FLXLP_inboard_spline(zz_idx), 'ro')
size(zz)
zz(zz_idx)



f_Bz_EF_z200 = fopen('180515009_9650_psi_r_z0.2.dat', 'r');
f_Bz_EF_z100 = fopen('180515009_9650_psi_r_z0.1.dat', 'r');
f_Bz_EF_z0 = fopen('180515009_9650_psi_r_z0.0.dat', 'r');
f_Bz_EF_z_nega_100 = fopen('180515009_9650_psi_r_z-0.1.dat', 'r');
f_Bz_EF_z_nega_200 = fopen('180515009_9650_psi_r_z-0.2.dat', 'r');
for i = 1 : 4
    temp = textscan(f_Bz_EF_z200,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1})
end
R_2d = zeros(1, 50)
psi = zeros(1, 50)
for i = 1 : 50
    temp = textscan(f_Bz_EF_z200,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1))
    psi(1, i) = str2double(temp_m(2))
end
plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(5));
Psi_z200 =psi + plus_inboard_flux;

for i = 1 : 4
    temp = textscan(f_Bz_EF_z100,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1})
end
R_2d = zeros(1, 50)
psi = zeros(1, 50)
for i = 1 : 50
    temp = textscan(f_Bz_EF_z100,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1))
    psi(1, i) = str2double(temp_m(2))
end
plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(4));
Psi_z100 =psi + plus_inboard_flux;

for i = 1 : 4
    temp = textscan(f_Bz_EF_z0,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
end
R_2d = zeros(1, 50);
psi = zeros(1, 50);
for i = 1 : 50
    temp = textscan(f_Bz_EF_z0,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1));
    psi(1, i) = str2double(temp_m(2));
end
plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(3));
Psi_z0 =psi + plus_inboard_flux;

for i = 1 : 4
    temp = textscan(f_Bz_EF_z_nega_100,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
end
R_2d = zeros(1, 50);
psi = zeros(1, 50);
for i = 1 : 50
    temp = textscan(f_Bz_EF_z_nega_100,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1));
    psi(1, i) = str2double(temp_m(2));
end
plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(2));
Psi_z_nega_100 =psi + plus_inboard_flux;

for i = 1 : 4
    temp = textscan(f_Bz_EF_z_nega_200,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
end
R_2d = zeros(1, 50);
psi = zeros(1, 50);
for i = 1 : 50
    temp = textscan(f_Bz_EF_z_nega_200,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    R_2d(1, i) = str2double(temp_m(1));
    psi(1, i) = str2double(temp_m(2));
end
plus_inboard_flux = ones(1, 50) * FLXLP_inboard_spline(zz_idx(1));
Psi_z_nega_200 =psi + plus_inboard_flux;





% figure('Name','Flux by 2D','NumberTitle','off')
% plot(R, Psi_z200)
% hold on
% plot(R, Psi_z100)
% hold on
% plot(R, Psi_z0)
% hold on
% plot(R, Psi_z_nega_100)
% hold on
% plot(R, Psi_z_nega_200)
% hold on
% plot(0.11, FLXLP_inboard_spline(zz_idx(5)), 'o')



% filename = '180515009_9650_psi_r_z0.2.dat';
% M = csvread(filename)


% psi_r = textscan(fid61,'%s',1,'Delimiter','\n');
% psi_r = strsplit(psi_r);
% c = textscan(f_Bz_EF, '%d %d'); % 4��̃f�[�^��int32�^�œǂݍ���
% pri_r = str2double(psi_r);
% disp(c)



if (IUTST <= 9)
    NFLX = 2*II;
elseif (IUTST == 10)
    NFLX = II;
elseif (IUTST == 11)
    NFLX = II;
end
fprintf('%s %d\r\n','NFLX=',NFLX);
%      
% ******************************************************
%   T-Probe
% ******************************************************
%
% frewind(fid61);% 
% textscan(fid61,'%s',1,'Delimiter','\n');  
frewind(fid68);% 
textscan(fid68,'%s',1,'Delimiter','\n');  
FACTR = 1.5D0;
II = 0;
%%C#####      LGB=1
LGB = 1;%21
LGBH = floor(LGB/2)+1;
LG = 0;
%%C####      MGB=1
MGB = 1; % �����ƂɃf�[�^���Ԉ����Ă��邩�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I
MGBH = floor(MGB/2)+1;
MG = 0;
for I = 1:IMAX_b
    temp = textscan(fid68,'%s',1,'Delimiter','\n');
    temp_m = strsplit(temp{1}{1});
    R = str2double(temp_m(1));
    Z = str2double(temp_m(2));
    %PSI = str2double(temp_m(3));
%     BZ = str2double(temp_m(4));
    BR = str2double(temp_m(5)); %      READ(61,*),R,Z,PSI,BZ,BR
    
    BZ = BZ_Distribution(1, I);

    
    %
    E = abs((R-RSENS)/RSENS);
    LG = LG + 1;
    LLG = LG;
    if (LG == LGB)  % 1�����Ń��Z�b�g����
        LG = 0;
    end
    MG = MG + 1;
    MMG = MG;
    if (MG == MGB) 
        MG = 0;
    end
    %
    %%C      IF(LLG.NE.LGBH) GOTO 200  
    if and (E < 1.0D-5, MMG ~= MGBH)
        continue  
    elseif and (E >= 1.0D-5, LLG ~= LGBH)
        continue

       elseif and(LWALL > 0, and (R < 0.12, abs(Z) < 0.9)); % �ǉ����̃Z���T�[����2016.9.5ushiki
           continue
    end      

    RR = R;
    ZZ = -Z;
    REND = RR - BR*FACTR;
    ZEND = ZZ + BZ*FACTR;
    fprintf(fid63,'%d %d\r\n',REND,ZEND);
    fprintf(fid63,'%d %d\r\n',RR,ZZ);
    fprintf(fid63,'%s\r\n','==');
%   
    II=II+1;
    if (IUTST <= 9)
        RS(NFLX+2*II-1) = R;
        RS(NFLX+2*II) = R;
        ZS(NFLX+2*II-1) = Z;     
        ZS(NFLX+2*II) = -Z;     
        Z1 = Z; BR1 = BR;
        Z2 = -Z; BR2 = -BR;
        TET1 = atan2(BZ,BR1);
        TET2 = atan2(BZ,BR2);
        % ���͐��̐ڐ������̎�����E���I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I            
        TET(NFLX+2*II-1) = TET1;%  !  
        TET(NFLX+2*II) = TET2;  %  !�@
        ITYPE(NFLX+2*II-1) = 1;
        ITYPE(NFLX+2*II) = 1;            
        BB1 = BR1*cos(TET1) + BZ*sin(TET1);
        BB2 = BR2*cos(TET2) + BZ*sin(TET2);
        XBR1 = BB1*cos(TET1);
        XBZ1 = BB1*sin(TET1);
        XBR2 = BB2*cos(TET2);
        XBZ2 = BB2*sin(TET2);
        fprintf(fid60,'%d %d %d %d %d %d %d %d\r\n',R,Z1,BB1,TET1,XBR1,BR1,XBZ1,BZ);
        EPS1R = abs((XBR1-BR1)/BR1);
        EPS1Z = abs((XBZ1-BZ)/BZ);
        fprintf(fid60,'%d %d %d %d %d %d %d %d\r\n',R,Z2,BB2,TET2,XBR2,BR2,XBZ2,BZ);
        EPS2R = abs((XBR2-BR2)/BR2);
        EPS2Z = abs((XBZ2-BZ)/BZ);
        if (EPS1R > 1.0E-7)
            fprintf(WAHAHA, '%d EPS1R= %d\r\n', I,EPS1R);
        end
        if (EPS1Z > 1.0E-7)
            fprintf(WAHAHA, '%d EPS1Z= %d\r\n', I,EPS1Z);   
        end
        if (EPS2R > 1.0E-7)
            fprintf(WAHAHA, '%d EPS2R= %d\r\n', I,EPS2R);
        end
        if (EPS2Z > 1.0E-7)
            fprintf(WAHAHA, '%d EPS2Z= %d\r\n', I,EPS2Z);
        end
        TPRB(2*II-1) = BB1;
        TPRB(2*II) = BB2;

    elseif (IUTST == 10)
        RS(NFLX+II) = R;
        ZS(NFLX+II) = Z;          
        Z1 = Z; BR1 = BR;
%         TET1 = atan2(BZ,BR1);       
         TET1 = 0.5 * pi;
        TET(NFLX+II) = TET1;%�@
        ITYPE(NFLX+II) = 1;            
        %����
%         BB1 = BR1*cos(TET1) + BZ*sin(TET1);
%         disp(BB1)
        BB1 = BZ;
        XBR1 = BB1*cos(TET1);
        XBZ1 = BB1*sin(TET1);
        fprintf(fid60,'%d %d %d %d %d %d %d %d\r\n',R,Z1,BB1,TET1,XBR1,BR1,XBZ1,BZ);
        EPS1R = abs((XBR1-BR1)/BR1);
        EPS1Z = abs((XBZ1-BZ)/BZ);
        if (EPS1R > 1.0E-7)
            fprintf(WAHAHA, '%d EPS1R= %d\r\n', I,EPS1R);
        end
        if (EPS1Z > 1.0E-7)
            fprintf(WAHAHA, '%d EPS1Z= %d\r\n', I,EPS1Z);   
        end
        TPRB(II) = BB1;
        elseif (IUTST == 11)
        RS(NFLX+II) = R;
        ZS(NFLX+II) = Z;          
        Z1 = Z; BR1 = BR;
        TET1 = atan2(BZ,BR1);
        % ���͐��̐ڐ������̎�����E���I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I            
        TET(NFLX+II) = TET1;%�@
        ITYPE(NFLX+II) = 1;            
        BB1 = BR1*cos(TET1) + BZ*sin(TET1);
        XBR1 = BB1*cos(TET1);
        XBZ1 = BB1*sin(TET1);
        fprintf(fid60,'%d %d %d %d %d %d %d %d\r\n',R,Z1,BB1,TET1,XBR1,BR1,XBZ1,BZ);
        EPS1R = abs((XBR1-BR1)/BR1);
        EPS1Z = abs((XBZ1-BZ)/BZ);
    end
end %200
if (IUTST <= 9)
    NTPB = 2*II;
elseif (IUTST == 10)
    NTPB = II;
elseif (IUTST == 11)
    NTPB = II;
end
fprintf('%s %d\r\n','NTPB=',NTPB);
%
% N-Probe	
%     Not defined.      
NNPB=0;
%%      NMAX=NFLX+NAPB
NMAX=NFLX+NTPB;
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++
fprintf(WAHAHA,'%s\r\n\r\n','Sensor data has been successfully installed!');
%++++++++++++++++++++++++++++++++++++++++++++++++++++
for I = 1:NMAX
	fprintf(WAHAHA,'%d Type=%d RS/ZS= %d %d Tet= %d\r\n',I,ITYPE(I),RS(I),ZS(I),TET(I));
end
%
fid110 = fopen('SENPOS0.txt','w'); % 110
fid111 = fopen('SENPOS1.txt','w'); % 111
fid112 = fopen('SENPOS2.txt','w'); % 112
for I = 1:NMAX
    if (ITYPE(I) == 0) 
        fprintf(fid110,'%d %d\r\n',RS(I),ZS(I));
    elseif (ITYPE(I) == 1)
        fprintf(fid111,'%d %d\r\n',RS(I),ZS(I));
    elseif(ITYPE(I) == 2)
        fprintf(fid112,'%d %d\r\n',RS(I),ZS(I));
    end
end
fclose(fid61);
fclose(fid60);
fclose(fid62);
fclose(fid63);
fclose(fid64);
fclose(fid65);
fclose(fid110);
fclose(fid111);
fclose(fid112);
%
%  *************************************************************************
%     Generation of CCS input data
%  *************************************************************************
%
fprintf('Generation of CCS input data ** START ***\n');
fid10 = fopen('CCSinput_UTST(temp).txt','w');
%           
fprintf(fid10,'%s\r\n','*');
fprintf(fid10,'%s\r\n','*** CCS Test input for UTST generated in PreUTST ***'); 
fprintf(fid10,'%s\r\n','** NTPB/NNPB/NFLX=(No. of T-/N-probes/Flux-Loops) **');
fprintf(fid10,'   %d     %d     %d\r\n',NTPB,NNPB,NFLX); 
fprintf(fid10,'%s\r\n','**  GETA (SSURF)');
fprintf(fid10,'  %s\r\n','0.0000E+00');
fprintf(fid10,'%s\r\n','* T-Probe');
for II = 1:NTPB
    fprintf(fid10,' %d\r\n',TPRB(II));
end
fprintf(fid10,'%s\r\n','* N-Probe');
fprintf(fid10,'%s\r\n','* Flux-Loop');
for II = 1:NFLX
    fprintf(fid10,' %d\r\n',FLXLP(II));
end
fprintf(fid10,'%s\r\n','****** MINR * MAXR * MINZ * MAXZ ****');
fprintf(fid10,'  %s\r\n','10   90  -100  100');
fprintf(fid10,'%s\r\n','*********');
fprintf(fid10,'%s\r\n','* ---�R�C���d���f�[�^�̕���---�P��[kA]');
fprintf(fid10,'%s\r\n','* EF');
fprintf(fid10,'%s\r\n','* PF#1');
fprintf(fid10,'%s\r\n','* PF#2');
fprintf(fid10,'%s\r\n','* PF#3');
fprintf(fid10,'%s\r\n','* PF#4');

PF_Current_for_CCS = zeros(4, 1);
for i = 1 : 4
    % �㉺�̃R�C���̓d���l�̕��ς�p����
    PF_Current_for_CCS(i, 1) = (Coil_Current(2 * i - 1, round(time_CCS / (0.5 * 0.001))) + Coil_Current(2 * i, round(time_CCS / (0.5 * 0.001)))) / 2;
    subplot(4, 1, i); plot(time, Coil_Current(2 * i - 1, :), time_CCS, Coil_Current(2 * i - 1, round(time_CCS / (0.5 * 0.001))), 'o', 'DisplayName', 'PF%d (upper side)')
    hold on
    subplot(4, 1, i); plot(time, Coil_Current(2 * i, :), 'DisplayName', 'PF%d (lower side)')
    legend
end
fprintf('%d', PF_Current_for_CCS)


% ���S�X�L�[�f�[�^���C���v�b�g�iEF�͎�ł��j
% fprintf(fid10,'%s\r\n', '-0.023');
% fprintf(fid10,'%s\r\n', '-0.023');
% fprintf(fid10,'%s\r\n', '-0.023');
% fprintf(fid10,'%s\r\n', '-0.023');
% fprintf(fid10,'%s\r\n', '-0.023');
% fprintf(fid10,'%s\r\n', '-0.023');
% fprintf(fid10,'%s\r\n', PF_Current_for_CCS(1, 1));
% fprintf(fid10,'%s\r\n', PF_Current_for_CCS(2, 1));
% fprintf(fid10,'%s\r\n', PF_Current_for_CCS(3, 1));
% fprintf(fid10,'%s\r\n', PF_Current_for_CCS(4, 1));

I_EF = -(0.849 * (1.19 * EF_voltage - 5.32) - 5.56);
% EF��200turn�R�C����6����ɂȂ��Ă���

fprintf(fid10,'%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10,'%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10,'%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10,'%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10,'%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10,'%d\r\n', 0.001 * I_EF / 6);
fprintf(fid10,'%d\r\n', PF_Current_for_CCS(1, 1));
fprintf(fid10,'%d\r\n', PF_Current_for_CCS(2, 1));
fprintf(fid10,'%d\r\n', PF_Current_for_CCS(3, 1));
fprintf(fid10,'%d\r\n', PF_Current_for_CCS(4, 1));

% fprintf(fid10,'%d\r\n', 0);
% fprintf(fid10,'%d\r\n', 0);
% fprintf(fid10,'%d\r\n', 0);
% fprintf(fid10,'%d\r\n', 0);
% fprintf(fid10,'%d\r\n', 0);
% fprintf(fid10,'%d\r\n', 0);
% fprintf(fid10,'%d\r\n', 0);
% fprintf(fid10,'%d\r\n', PF_Current_for_CCS(2, 1));
% fprintf(fid10,'%d\r\n', 0);
% fprintf(fid10,'%d\r\n', 0);



if and(IUTST <= 10, IUTST ~= 8)
% EF���\���l��
%�E��n
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','2.40');
% fprintf(fid10,'   %s\r\n','-23.0');
% fprintf(fid10,'   %s\r\n','2.33');
% fprintf(fid10,'   %s\r\n','-22.2');

% 9650us
%     fprintf(fid10,'   %s\r\n','-0.28');
%     fprintf(fid10,'   %s\r\n','2.40');
%     fprintf(fid10,'   %s\r\n','-23.0');
%     fprintf(fid10,'   %s\r\n','2.33');
%     fprintf(fid10,'   %s\r\n','-22.2');
%%9500us
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','-13.0');
% fprintf(fid10,'   %s\r\n','28.4');
% fprintf(fid10,'   %s\r\n','-2.74');
% fprintf(fid10,'   %s\r\n','27.9');
%%9530us
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-2.3');
% fprintf(fid10,'   %s\r\n','-27.2');
% fprintf(fid10,'   %s\r\n','3.93');
% fprintf(fid10,'   %s\r\n','-24.1');
%9450us
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','-0.047');
% fprintf(fid10,'   %s\r\n','1.37');
% fprintf(fid10,'   %s\r\n','-25.5');
% fprintf(fid10,'   %s\r\n','2.02');
% fprintf(fid10,'   %s\r\n','-23.2');
%%8000us
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0.047');
% fprintf(fid10,'   %s\r\n','0');
% fprintf(fid10,'   %s\r\n','0');
% fprintf(fid10,'   %s\r\n','2.64');
% fprintf(fid10,'   %s\r\n','0');
elseif (IUTST == 8)
    COIL_C = load(strcat('UTST_8/coil_eddy_',num2str(eddy_time),'.txt'));
    fprintf(fid10,'   %d\r\n',0.0);
    fprintf(fid10,'   %d\r\n',-COIL_C(1));
    fprintf(fid10,'   %d\r\n',-COIL_C(2));
    fprintf(fid10,'   %d\r\n',-COIL_C(3));
    fprintf(fid10,'   %d\r\n',-COIL_C(4));
end
fclose(fid10);
fprintf('Generation of CCS input data ** END ***\n');
end
%
%% INPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [CR,CZ,FFDAT,R0,Z0,RR,NCOUNT,RC,ZC,ECI,NTPB,NNPB,NAPB,NFLX,FLXLP,...
BSNSR,CAPPER,REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,IDECCS,SOU,REVN,ZEVN,GETA,NE,...
KCMX,ISHT,NINT,R_plot,Z_plot,CCS,MSEC]...
= INPUT(n100,Nedp,IUTST,MXINT,WAHAHA,NONC,AUTO,AUTOINPUT,xu,yu,xl,yl)
CR = zeros(1,MXINT); 
CZ = zeros(1,MXINT);  
FFDAT = zeros(1,n100); 
FLXLP = zeros(1,n100); 
BSNSR = zeros(1,n100);
REV = zeros(1,Nedp); 
ZEV = zeros(1,Nedp); 
RES = zeros(1,Nedp); 
ZES = zeros(1,Nedp); 
REVN = zeros(1,Nedp); 
ZEVN = zeros(1,Nedp);
R0 = 0;
Z0 = 0;
RR = 0;
CAPPER = 0;
IDECCS = 0;
SOU = 0;

NE = 0;
MSEC = 0; %ushiki
ISHT = 0;
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX   
if (IUTST > 0)
    INP = fopen('CCSinput_UTST(temp).txt','r');     
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      
else
    fprintf('�v���Y�}�̃^�C�v�F�@ Parabolic (1)�H\n');
    fprintf('                   Hollow    (2)�H\n');
    fprintf('                   Peaked    (3)�H\n');
    prompt = '                   Broad     (4)�H\n';
    KRTYPE = input(prompt);
%  Plasma Current Type: (KRTYPE=(1:parabolic, 2:hollow, 3:peaked 4:broad)) 
    if (KRTYPE == 1)
        INP = fopen('CCSinput_Parabolic.txt','r');
        fprintf('***** Parabolic Type has been chosen.\n');
    elseif (KRTYPE == 2)
        INP = fopen('CCSinput_Hollow.txt','r');
        fprintf('***** Hollow Type has been chosen.\n');
    elseif (KRTYPE == 3)
        INP = fopen('CCSinput_Peaked.txt','r');
        fprintf('***** Peaked Type has been chosen.\n');
    elseif (KRTYPE == 4)
        INP = fopen('CCSinput_Broad.txt','r');
        fprintf('***** Broad Type has been chosen.\n');
    end
end
HEAD = textscan(INP,'%1s',1); %801
while (HEAD{1}{1} == '*')
    HEAD = textscan(INP,'%1s',1); %801
end
temp = textscan(INP,'%s',1,'Delimiter','\n');
TITLE = strcat(HEAD{1}{1},temp{1}{1});
fprintf(WAHAHA, 'Project Name ==>   %s\r\n\r\n\r\n\r\n',TITLE);%      WRITE(IPR,250) TITLE
HEAD = textscan(INP,'%1s',1); 
while (HEAD{1}{1} == '*');
    textscan(INP,'%s',1,'Delimiter','\n'); %802
    HEAD = textscan(INP,'%1s',1); 
end
temp = textscan(INP,'%s',1,'Delimiter','\n');
Sensnum = strcat(HEAD{1}{1},temp{1}{1});
temp_m = strsplit(Sensnum);
NTPB = str2double(temp_m(1));
NNPB = str2double(temp_m(2));
NFLX = str2double(temp_m(3));
NAPB = NTPB + NNPB;
%
HEAD = textscan(INP,'%1s',1); 
while (HEAD{1}{1} == '*')
    textscan(INP,'%s',1,'Delimiter','\n'); %803
    HEAD = textscan(INP,'%1s',1); 
end
temp = textscan(INP,'%s',1,'Delimiter','\n');
GETA = strcat(HEAD{1}{1},temp{1}{1});
%
for I = 1:NAPB+NFLX
    HEAD = textscan(INP,'%1s',1); 
    while (HEAD{1}{1} == '*')
        textscan(INP,'%s',1,'Delimiter','\n'); %804
        HEAD = textscan(INP,'%1s',1); 
    end
    temp = textscan(INP,'%s',1,'Delimiter','\n');
    FFDAT(I) = str2double(strcat(HEAD{1}{1},temp{1}{1}));
    if (I <= NAPB)
    	FFDAT(I)=FFDAT(I);  
    elseif (IUTST > 0)
	    FFDAT(I) = FFDAT(I)/(2*pi);%      ! (����) �����Ǝ����֐��̈Ⴂ���l��
    end
%C      FFDAT = MEASURED VALUE // T-Probe, N-Probe, FLUX LOOP
end
%*****
if (AUTO == 0)
    fprintf('Input the standard deviation of Gaussian noise\n');
    prompt = '(e.g.,   0.03 for 3% sigma; 0.00 for no noise)\n';
    SIGM = input(prompt); 
    if (SIGM ~= 0)
        prompt = 'Choose the SEED number of Gaussian noise (1,2,3...)\n';
        SEED = input(prompt);
    else
        SEED = 1;
    end
else
    SIGM = AUTOINPUT(4);
    SEED = AUTOINPUT(5);
end
%IDUM = -100; %-100
rng(SEED); % ���񓯂������𔭐�������I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I
%rng(2); % ���񓯂������𔭐�������I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I
%rng('shuffle'); % �V���b�t���I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I
GASDEV = randn(NAPB+NFLX,1); % ���K����
%plot(GASDEV)

%�K�E�V�A���m�C�Y��t��
for I=1:NAPB+NFLX
    FFDAT(I) = FFDAT(I)*(1.0+SIGM*GASDEV(I,1));
end
for I=1:NAPB
	BSNSR(I) = FFDAT(I);
end
%
for I=1:NFLX;
	FLXLP(I) = FFDAT(I+NAPB);
end
%
HEAD = textscan(INP,'%1s',1); 
while (HEAD{1}{1} == '*')
    textscan(INP,'%s',1,'Delimiter','\n'); %805
    HEAD = textscan(INP,'%1s',1); 
end
temp = textscan(INP,'%s',1,'Delimiter','\n');
temp_m = strcat(HEAD{1}{1},temp{1}{1});
temp_1 = strsplit(temp_m);
MINR = str2double(temp_1(1));
MAXR = str2double(temp_1(2));
MINZ = str2double(temp_1(3));
MAXZ = str2double(temp_1(4)); %      READ(61,*),R,Z,PSI,BZ,BR
fprintf('    %d     %d     %d     %d\r\n',MINR,MAXR,MINZ,MAXZ);
%   
%** CR(),CZ() ==> Data points to make the contour map of magnetic flux  
ICRE = 5.*(IUTST == 0) + 1.*(IUTST ~= 0);
JCRE = 5.*(IUTST == 0) + 2.*(IUTST ~= 0);
NINT = 0;
for I = MINR:ICRE:MAXR
      NCOUNT = 0;
      CCR = I/100.0;
      for J = MINZ:JCRE:MAXZ
          NINT = NINT+1;
          NCOUNT = NCOUNT+1;
          CR(NINT) = CCR;
          CZ(NINT) = J/100.0;
      end
end
R_plot = round((MAXR -MINR)/ICRE);
Z_plot = round((MAXZ -MINZ)/JCRE);
%plot(CR,CZ,'o');
fprintf('%s %d\r\n','NINT =',NINT);
%%%
if (AUTO == 0)
    prompt = 'CCS�̐���?\n';
    CCS = input(prompt);
else
    CCS = AUTOINPUT(6);
end
if (CCS > 0)
    if (AUTO == 0)
        prompt = 'CCS�`��͑ȉ~�HD�^�H=(0/1)\n';
        IDECCS = input(prompt);
    else
        IDECCS = AUTOINPUT(7);
    end
    if (CCS == 1)
        if (AUTO == 0)
            fprintf('CCS�̒��S(R,Z)�̏����l��?  (��F (0.35, 0.0))\n');
            prompt = 'R = ';
            R0 = input(prompt);
            prompt = 'Z = ';
            Z0 = input(prompt);
            prompt = 'CCS��R�����̔��a��?   (��F 0.04)\n';
            RR = input(prompt);
           %%
            %! �R�[�V�[�����ʂ̔�~�`�x��
            prompt = 'CCS�̏c����(�c�a��Z/���a��R)��? (��F 1.3)\n';
            CAPPER = input(prompt);
           %%
            if (IDECCS == 0)
                prompt = '��K�����E�v�f��(CC�ߓ_����1/3)��?�i�����l��6�j\n';
                NE = input(prompt);
            else
                fprintf('D�^CCS�̋��E�v�f��(CC�ߓ_����1/3)��D_CCS �ŗ^���܂�\n');
                prompt = '��K�����E�v�f��(CC�ߓ_����1/3)D_CCS�ȗ������?�i�����l��2�j\n';
                MSEC(1) = input(prompt);
                prompt = '��K�����E�v�f��(CC�ߓ_����1/3)D_CCS�ȗ�������?�i�����l��2�j\n';
                MSEC(2) = input(prompt);
                prompt = '��K�����E�v�f��(CC�ߓ_����1/3)D_CCS�C���{�[�h����?�i�����l��2�j\n';
                MSEC(3) = input(prompt);
            end
        %%
            prompt = '##\r\nCCS�`��C���̍ہALCMS��CCS�̑�����́H�i��F=5.0�j\n';
            SOU = input(prompt);
        else
            R0 =     xu;
            Z0 =     yu;
            RR =     AUTOINPUT(10);
            CAPPER = AUTOINPUT(11);
            NE =     AUTOINPUT(12);
            MSEC(1) =AUTOINPUT(13);
            MSEC(2) =AUTOINPUT(14);
            MSEC(3) =AUTOINPUT(15);            
            SOU =    AUTOINPUT(16);
        end
    elseif (CCS > 1)
        R0 = zeros(1,CCS);
        Z0 = zeros(1,CCS);
        RR = zeros(1,CCS);
        CAPPER = zeros(1,CCS);
        NE = zeros(1,CCS);
        SOU = zeros(1,CCS);
    %     prompt = 'CCS�`��͑ȉ~�HD�^�H=(0/1)\n';
    %     IDECCS = input(prompt);
        for I = 1:CCS
            if (AUTO == 0)
                fprintf('%d%s\n',I,'�Ԗڂ�CCS�̒��S(R,Z)�̏����l��?  (��F (0.35, 0.0))');
                prompt = 'R = ';
                R0(I) = input(prompt);
                prompt = 'Z = ';
                Z0(I) = input(prompt);
                prompt = strcat(num2str(I),'�Ԗڂ�CCS��R�����̔��a��?   (��F 0.04)\n');
                RR(I) = input(prompt);
               %%
                %! �R�[�V�[�����ʂ̔�~�`�x��
                prompt = strcat(num2str(I),'�Ԗڂ�CCS�̏c����(�c�a��Z/���a��R)��? (��F 1.3)\n');
                CAPPER(I) = input(prompt);
               %%
                if (IDECCS == 0)
                    prompt = strcat(num2str(I),'�Ԗڂ�CCS��K�����E�v�f��(CC�ߓ_����1/3)��?�i�����l��6�j\n');
                    NE(I) = input(prompt);
                else
                    fprintf('D�^CCS�̋��E�v�f��(CC�ߓ_����1/3)��D_CCS �ŗ^���܂�\n');
                    prompt = strcat(num2str(I),'��K�����E�v�f��(CC�ߓ_����1/3)D_CCS�ȗ������?�i�����l��2�j\n');
                    MSEC(I,1) = input(prompt);
                    prompt = strcat(num2str(I),'��K�����E�v�f��(CC�ߓ_����1/3)D_CCS�ȗ�������?�i�����l��2�j\n');
                    MSEC(I,2) = input(prompt);
                    prompt = strcat(num2str(I),'��K�����E�v�f��(CC�ߓ_����1/3)D_CCS�C���{�[�h����?�i�����l��2�j\n');
                    MSEC(I,3) = input(prompt);
                end
    %%
                prompt = strcat(num2str(I),'##\r\nCCS�`��C���̍ہALCMS��CCS�̑�����́H�i��F=5.0�j\n');
                SOU(I) = input(prompt);
            else
                if I == 1
                    R0(I) =     xl;
                    Z0(I) =     yl;
                elseif I == 2
                    R0(I) =     xu;
                    Z0(I) =     yu;                    
                end
                RR(I) =     AUTOINPUT(10 + (I-1)*9);
                CAPPER(I) = AUTOINPUT(11 + (I-1)*9);
                NE(I) =     AUTOINPUT(12 + (I-1)*9);
                MSEC(I,1) = AUTOINPUT(13 + (I-1)*9);
                MSEC(I,2) = AUTOINPUT(14 + (I-1)*9);
                MSEC(I,3) = AUTOINPUT(15 + (I-1)*9);            
                SOU(I) =    AUTOINPUT(16 + (I-1)*9);                
            end
        end
    end
end
%
% �R�C���ʒu�f�[�^�̓ǂݎ��
if (IUTST > 0)
%     [KCMX,RC,ZC,ECI] = COIL_UTST(INP,WAHAHA); %OK
    [KCMX,RC,ZC,ECI] = COIL_UTST_2(INP,WAHAHA); %OK
else      
    ISHT=3;
%    CALL COIL60K(ISHT,INP,IPR,KCMX,RC,ZC,ECI,n900)
end
% ���C�Z���T�[�ʒu�f�[�^�̓ǂݎ��      
if (IUTST <= 0)          
%      CALL SENS60K(IPR,NMAX,RS,ZS,ITYPE,TET,nsemx)
end
%
% ���̏�̉Q�d���ߓ_�ʒu�f�[�^�̓ǂݎ�� 
if (AUTO == 0)
    prompt = '### ���̏�̉Q�d�����l������H (Yes/No)=(1/0)\n';   
    LCOND = input(prompt);
else
    LCOND = AUTOINPUT(26);
end
KNE=0;
KSE=0;
KNN=0;
KSN=0;
if (LCOND > 0) 
    if (IUTST > 0)
        [KNE,KNN,REV,ZEV,RES,ZES,KSE,KSN,REVN,ZEVN] = CONDSHL_UTST(WAHAHA,Nedp,NONC); %OK
%    else
%        CALL CONDSHL(REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,Nedp)
    end
%else
end
end
%
%% CONDSHL_UTST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ���̏�̉Q�d���ߓ_�ʒu�f�[�^�̐ݒ�Ɠǂݎ��
% (REVN, ZEVN)�F�Q�d���m�[�h���W
function [KNE,KNN,REV,ZEV,RES,ZES,KSE,KSN,REVN,ZEVN] = CONDSHL_UTST(WAHAHA,Nedp,NONC)
eddyNode = 7;% 3 standard %7mid�t��outboard�ɒu���Ȃ�
if (eddyNode == 1)
    MAXM = 13;
    Z = 0.15;
    %Z=0.1; Q=0.03;
    %Z=0.050 ; Q=0.03;%   ! ����ʖځA�L����
    %Z=0.030; Q=0.03;   % Very standard!!
    %Z=0.001 ; Q=0.001;%   ! ����ʖځA������
    %C Z=0.010 ; Q=0.03;   ! ����ʖ�
    REV = zeros(1,Nedp); 
    ZEV = zeros(1,Nedp); 
    RES = zeros(1,Nedp); 
    ZES = zeros(1,Nedp);
    RSEC = zeros(1,MAXM+1); 
    ZSEC = zeros(1,MAXM+1); 
    REVN = zeros(1,Nedp); 
    ZEVN = zeros(1,Nedp);
    INBOAD = 0.5;
    %       1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    MSEC = [1, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 1];
    RSEC(1)= 0.694;        ZSEC(1)= 0.0;
    RSEC(2)= 0.694;        ZSEC(2)= 0.285;
    RSEC(3)= 0.5985;       ZSEC(3)= 0.285;
    RSEC(4)= 0.5985;       ZSEC(4)= 0.285+Z;
    RSEC(5)= 0.5985;       ZSEC(5)= 0.9985;
    RSEC(6)= 0.10815;      ZSEC(6)= 0.9985;
    RSEC(7)= 0.10815;      ZSEC(7)= INBOAD;
    RSEC(8)= 0.10815;      ZSEC(8)= -INBOAD;
    RSEC(9)= 0.10815;      ZSEC(9)= -0.9985;
    RSEC(10)= 0.5985;      ZSEC(10)= -0.9985;
    RSEC(11)=  0.5985;     ZSEC(11)= -0.285-Z;
    RSEC(12)= 0.5985;      ZSEC(12)= -0.285;
    RSEC(13)= 0.694;       ZSEC(13)= -0.285;
    RSEC(14)= 0.694;       ZSEC(14)= 0.0;   
elseif(eddyNode == 2)
    MAXM = 13;
    Z = 0.15;
    REV = zeros(1,Nedp); 
    ZEV = zeros(1,Nedp); 
    RES = zeros(1,Nedp); 
    ZES = zeros(1,Nedp);
    RSEC = zeros(1,MAXM+1); 
    ZSEC = zeros(1,MAXM+1); 
    REVN = zeros(1,Nedp); 
    ZEVN = zeros(1,Nedp);
    INBOAD = 0.5;
    %       1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    MSEC = [1, 1, 1, 3, 1, 2, 5, 2, 1, 3, 1, 1, 1];
    RSEC(1)= 0.694;        ZSEC(1)= 0.0;
    RSEC(2)= 0.694;        ZSEC(2)= 0.285;
    RSEC(3)= 0.5985;       ZSEC(3)= 0.285;
    RSEC(4)= 0.5985;       ZSEC(4)= 0.285+Z;
    RSEC(5)= 0.5985;       ZSEC(5)= 0.9985;
    RSEC(6)= 0.10815;      ZSEC(6)= 0.9985;
    RSEC(7)= 0.10815;      ZSEC(7)= INBOAD;
    RSEC(8)= 0.10815;      ZSEC(8)= -INBOAD;
    RSEC(9)= 0.10815;      ZSEC(9)= -0.9985;
    RSEC(10)= 0.5985;      ZSEC(10)= -0.9985;
    RSEC(11)=  0.5985;     ZSEC(11)= -0.285-Z;
    RSEC(12)= 0.5985;      ZSEC(12)= -0.285;
    RSEC(13)= 0.694;       ZSEC(13)= -0.285;
    RSEC(14)= 0.694;       ZSEC(14)= 0.0;   
elseif(eddyNode == 3)
    MAXM = 13;
    Z = 0.15;
    REV = zeros(1,Nedp); 
    ZEV = zeros(1,Nedp); 
    RES = zeros(1,Nedp); 
    ZES = zeros(1,Nedp);
    RSEC = zeros(1,MAXM+1); 
    ZSEC = zeros(1,MAXM+1); 
    REVN = zeros(1,Nedp); 
    ZEVN = zeros(1,Nedp);
    INBOAD = 0.5;
    %       1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    MSEC = [1, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 1];
    RSEC(1)= 0.694;        ZSEC(1)= 0.0;
    RSEC(2)= 0.694;        ZSEC(2)= 0.285;
    RSEC(3)= 0.5985;       ZSEC(3)= 0.285;
    RSEC(4)= 0.5985;       ZSEC(4)= 0.285+Z;
    RSEC(5)= 0.5985;       ZSEC(5)= 0.9985;
    RSEC(6)= 0.10815;      ZSEC(6)= 0.9985;
    RSEC(7)= 0.10815;      ZSEC(7)= INBOAD;
    RSEC(8)= 0.10815;      ZSEC(8)= -INBOAD;
    RSEC(9)= 0.10815;      ZSEC(9)= -0.9985;
    RSEC(10)= 0.5985;      ZSEC(10)= -0.9985;
    RSEC(11)=  0.5985;     ZSEC(11)= -0.285-Z;
    RSEC(12)= 0.5985;      ZSEC(12)= -0.285;
    RSEC(13)= 0.694;       ZSEC(13)= -0.285;
    RSEC(14)= 0.694;       ZSEC(14)= 0.0;   
    
%     plot(RSEC, ZSEC, 'x')
%     
    
    elseif (eddyNode == 4)
    MAXM = 13;
    Z = 0.15;
    REV = zeros(1,Nedp); 
    ZEV = zeros(1,Nedp); 
    RES = zeros(1,Nedp); 
    ZES = zeros(1,Nedp);
    RSEC = zeros(1,MAXM+1); 
    ZSEC = zeros(1,MAXM+1); 
    REVN = zeros(1,Nedp); 
    ZEVN = zeros(1,Nedp);
    INBOAD = 0.5;
    %       1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    MSEC = [1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1];
    RSEC(1)= 0.697;        ZSEC(1)= 0.0;
    RSEC(2)= 0.697;        ZSEC(2)= 0.2925;
    RSEC(3)= 0.59925;       ZSEC(3)= 0.2925;
    RSEC(4)= 0.59925;       ZSEC(4)= 0.2925+Z;
    RSEC(5)= 0.59925;       ZSEC(5)= 0.99925;
    RSEC(6)= 0.10615;      ZSEC(6)= 0.99925;
    RSEC(7)= 0.10615;      ZSEC(7)= INBOAD;
    RSEC(8)= 0.10615;      ZSEC(8)= -INBOAD;
    RSEC(9)= 0.10615;      ZSEC(9)= -0.99925;
    RSEC(10)= 0.59925;      ZSEC(10)= -0.99925;
    RSEC(11)=  0.59925;     ZSEC(11)= -0.2925-Z;
    RSEC(12)= 0.59925;      ZSEC(12)= -0.2925;
    RSEC(13)= 0.697;       ZSEC(13)= -0.2925;
    RSEC(14)= 0.697;       ZSEC(14)= 0.0;   
    elseif (eddyNode == 5)
    MAXM = 15;
    Z = 0.15;
    REV = zeros(1,Nedp); 
    ZEV = zeros(1,Nedp); 
    RES = zeros(1,Nedp); 
    ZES = zeros(1,Nedp);
    RSEC = zeros(1,MAXM+1); 
    ZSEC = zeros(1,MAXM+1); 
    REVN = zeros(1,Nedp); 
    ZEVN = zeros(1,Nedp);
    INBOAD = 0.1;
    %       1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    MSEC = [1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1];
    RSEC(1)= 0.697;        ZSEC(1)= 0.0;
    RSEC(2)= 0.697;        ZSEC(2)= 0.2925;
    RSEC(3)= 0.59925;       ZSEC(3)= 0.2925;
    RSEC(4)= 0.59925;       ZSEC(4)= 0.2925+Z;
    RSEC(5)= 0.59925;       ZSEC(5)= 0.99925;
    RSEC(6)= 0.130523;      ZSEC(6)= 0.99925;
    RSEC(7)= 0.10615;      ZSEC(7)= 0.99925;
    RSEC(8)= 0.10615;      ZSEC(8)= 0.99925-INBOAD;
    RSEC(9)= 0.10615;      ZSEC(9)= -0.99925+INBOAD;
    RSEC(10)= 0.10615;      ZSEC(10)= -0.99925;
    RSEC(11)= 0.130523;      ZSEC(11)= -0.99925;
    RSEC(12)= 0.59925;      ZSEC(12)= -0.99925;
    RSEC(13)=  0.59925;     ZSEC(13)= -0.2925-Z;
    RSEC(14)= 0.59925;      ZSEC(14)= -0.2925;
    RSEC(15)= 0.697;       ZSEC(15)= -0.2925;
    RSEC(16)= 0.697;       ZSEC(16)= 0.0;   
    elseif (eddyNode == 6)
    MAXM = 11;
    Z = 0.15;
    REV = zeros(1,Nedp); 
    ZEV = zeros(1,Nedp); 
    RES = zeros(1,Nedp); 
    ZES = zeros(1,Nedp);
    RSEC = zeros(1,MAXM+1); 
    ZSEC = zeros(1,MAXM+1); 
    REVN = zeros(1,Nedp); 
    ZEVN = zeros(1,Nedp);
    INBOAD = 0.5;
    %       1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
%     MSEC = [1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1];
    MSEC = [1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1];
    RSEC(1)= 0.697;        ZSEC(1)= 0.0;
    RSEC(2)= 0.697;        ZSEC(2)= 0.2925;
    RSEC(3)= 0.59925;       ZSEC(3)= 0.2925;
    RSEC(4)= 0.59925;       ZSEC(4)= 0.2925+Z;
    RSEC(5)= 0.59925;       ZSEC(5)= 0.99925;
    RSEC(6)= 0.10615;      ZSEC(6)= 0.99925;
    RSEC(7)= 0.10615;      ZSEC(7)= -0.99925;
    RSEC(8)= 0.59925;      ZSEC(8)= -0.99925;
    RSEC(9)=  0.59925;     ZSEC(9)= -0.2925-Z;
    RSEC(10)= 0.59925;      ZSEC(10)= -0.2925;
    RSEC(11)= 0.697;       ZSEC(11)= -0.2925;
    RSEC(12)= 0.697;       ZSEC(12)= 0.0;   
 elseif(eddyNode == 7)
    MAXM = 11;
    Z = 0.2;
    Z = 0.4;
    REV = zeros(1,Nedp); 
    ZEV = zeros(1,Nedp); 
    RES = zeros(1,Nedp); 
    ZES = zeros(1,Nedp);
    RSEC = zeros(1,MAXM+1); 
    ZSEC = zeros(1,MAXM+1); 
    REVN = zeros(1,Nedp); 
    ZEVN = zeros(1,Nedp);
    INBOAD = 0.5;
    %       1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
%     MSEC = [1, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 1];
%     MSEC = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    MSEC = [1, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1, 0, 0];
%     RSEC(1)= 0.694;        ZSEC(1)= 0.0;
    RSEC(1)= 0.694;        ZSEC(1)= 0.285;
    RSEC(2)= 0.5985;       ZSEC(2)= 0.285;
    RSEC(3)= 0.5985;       ZSEC(3)= 0.285+Z;
    RSEC(4)= 0.5985;       ZSEC(4)= 0.9985;
    RSEC(5)= 0.10815;      ZSEC(5)= 0.9985;
    RSEC(6)= 0.10815;      ZSEC(6)= INBOAD;
    RSEC(7)= 0.10815;      ZSEC(7)= -INBOAD;
    RSEC(8)= 0.10815;      ZSEC(8)= -0.9985;
    RSEC(9)= 0.5985;      ZSEC(9)= -0.9985;
    RSEC(10)=  0.5985;     ZSEC(10)= -0.285-Z;
    RSEC(11)= 0.5985;      ZSEC(11)= -0.285;
    RSEC(12)= 0.694;       ZSEC(12)= -0.285;
%     RSEC(13)= 0.694;       ZSEC(13)= 0.285;
%     RSEC(14)= 0.694;       ZSEC(14)= 0.285;   
 elseif(eddyNode == 8)
    MAXM = 13;
%     Z = 0.2;
    REV = zeros(1,Nedp); 
    ZEV = zeros(1,Nedp); 
    RES = zeros(1,Nedp); 
    ZES = zeros(1,Nedp);
    RSEC = zeros(1,MAXM+1); 
    ZSEC = zeros(1,MAXM+1); 
    REVN = zeros(1,Nedp); 
    ZEVN = zeros(1,Nedp);
    INBOAD = 0.5;
    %       1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
%     MSEC = [1, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 1];
%     MSEC = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    MSEC = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0];
%     RSEC(1)= 0.694;        ZSEC(1)= 0.0;
    RSEC(1)= 0.694;        ZSEC(1)= 0.285;
    RSEC(2)= 0.5985;       ZSEC(2)= 0.285;
    RSEC(3)= 0.5985;       ZSEC(3)= 0.285+(0.9985 - 0.285) / 3;
    RSEC(4)= 0.5985;       ZSEC(4)= 0.285+(0.9985 - 0.285) / 3 * 2;
    RSEC(5)= 0.5985;       ZSEC(5)= 0.9985;
    RSEC(6)= 0.10815;      ZSEC(6)= 0.9985;
    RSEC(7)= 0.10815;      ZSEC(7)= INBOAD;
    RSEC(8)= 0.10815;      ZSEC(8)= -INBOAD;
    RSEC(9)= 0.10815;      ZSEC(9)= -0.9985;
    RSEC(10)= 0.5985;      ZSEC(10)= -0.9985;
    RSEC(11)=  0.5985;     ZSEC(11)= -0.285-(0.9985 - 0.285) / 3 * 2;
    RSEC(12)=  0.5985;     ZSEC(12)= -0.285-(0.9985 - 0.285) / 3;
    RSEC(13)= 0.5985;      ZSEC(13)= -0.285;
    RSEC(14)= 0.694;       ZSEC(14)= -0.285;
%     RSEC(13)= 0.694;       ZSEC(13)= 0.285;
%     RSEC(14)= 0.694;       ZSEC(14)= 0.285;   
end
% �^��e���̉Q�d���ߓ_�̐ݒ�
% KNE=�^��e��̕������E�v�f��,  KNM=�^��e���̃��b�V���_��,  KNN=�^��e���̐ߓ_��
fid110 = fopen('VacuumVesselSegments.txt','w');
myformat = '%d %d\r\n';
for I = 1:MAXM+1
    fprintf(fid110,myformat,RSEC(I),ZSEC(I));
end
fclose(fid110);
II=1;
REV(II)=RSEC(1);
ZEV(II)=ZSEC(1);
KNE=0;
for I = 1:MAXM
    KNE = KNE+MSEC(I);
    CM = 2*MSEC(I);
    DELR = (RSEC(I+1)-RSEC(I))/CM;
    DELZ = (ZSEC(I+1)-ZSEC(I))/CM;
    for J = 1:2*MSEC(I)
        II = II+1;
        REV(II) = REV(II-1)+DELR;
        ZEV(II) = ZEV(II-1)+DELZ;
    end
end
%
KNM = II-1;
fid112 = fopen('VacuumVesselMeshPoints.txt','w');
for II = 1:KNM+1
       fprintf(fid110,myformat,REV(II),ZEV(II));
end  
%% �I�I�I�I�I�I�I�@ ��K���v�f�ߓ_�v�Z�̗\��n�@�@�@�I�I�I�I�I�I�I
% �^��e���̔�K���v�f�ߓ_���W�̍쐬
if (NONC > 0)
    KNN = KNE*3;
        I = 1:KNE;
	    REVN(3*I-2) = (5.D0*REV(2*I-1)+5.D0.*REV(2*I)-REV(2*I+1))/9.D0;
	    ZEVN(3*I-2) = (5.D0*ZEV(2*I-1)+5.D0.*ZEV(2*I)-ZEV(2*I+1))/9.D0;
	    REVN(3*I-1) = REV(2*I);
	    ZEVN(3*I-1) = ZEV(2*I);
	    REVN(3*I)   = (5.D0*REV(2*I+1)+5.D0.*REV(2*I)-REV(2*I-1))/9.D0;
	    ZEVN(3*I)   = (5.D0*ZEV(2*I+1)+5.D0.*ZEV(2*I)-ZEV(2*I-1))/9.D0;
%    
else
    KNN = KNM;
end
%
fprintf(WAHAHA,'%s\r\n','Position of Eddy current computing points');
fid110 = fopen('VacuumVesselNodePoints.txt','w'); %113
if (NONC == 0)
    for II = 1:KNM+1
        fprintf(WAHAHA,'%d %d\r\n',REV(II),ZEV(II));
        fprintf(fid110,'%d %d\r\n',REV(II),ZEV(II));
    end
        fprintf(WAHAHA,'KNE=%d    KNM=%d\r\n',KNE,KNM);
        fprintf('KNE=%d    KNM=%d\r\n',KNE,KNM);
        fprintf(WAHAHA,'%s\r\n','The end meshpoint agrees with the start meshpoint');
        fprintf('%s\r\n','The end meshpoint agrees with the start meshpoint');
else
    for II = 1:KNN
        fprintf(WAHAHA,'%d %d\r\n',REVN(II),ZEVN(II));
        fprintf(fid110,'%d %d\r\n',REVN(II),ZEVN(II));
    end
    fprintf(WAHAHA,'KNE/KNM/KNN = %d  %d  %d\r\n',KNE,KNM,KNN);
    fprintf('KNE/KNM/KNN = %d  %d  %d\r\n',KNE,KNM,KNN);
end 
%
%
% ���艻��̉Q�d���ߓ_�̐ݒ�      
% KSE=���艻�̕������E�v�f��,  KSN=���艻��̐ߓ_��
%***
%cd      KSE=1 
%prompt = '�d���V�[�g��ɋ��E�v�f��z�u����H (Yes/No)=(1/0)\n';   
%KSE = input(prompt);
KSE = 0; 
%***
KSN = KSE*2 + 1;
if (KSE == 0)
    KSN=0;
end
if (KSE > 0)
    fid112 = fopen('StabilizerPoints.txt','w'); % 114
    RES(1) = 0.25;
    RES(2) = 0.28;
    RES(3) = 0.31;
    ZES(1) = 0;     
    ZES(2) = 0;    
    ZES(3) = 0;     
    for I = 1:KSN
        fprintf(WAHAHA,'%d %d %d\r\n',I,RES(I),ZES(I));
        fprintf(fid112,'%d %d\r\n',RES(I),ZES(I));
    end
else
end 
fclose(fid110);
fclose(fid112);
end
%
%% FORM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% AA*X=FF
function [FC,BR,BZ,PSIFLX,PSIC,AMYU0,AA,FF,GETA,NCCS,fid99,fid100]= FORM(WAHAHA,...
IT,OLDGT,GETA,FF,AA,n100,n50,ECI,KCMX,RS,ZS,RC,ZC,ITYPE,TET,NTPB,NNPB,...
NAPB,NFLX,RCCN,ZCCN,NCCN,NCCS,RCCS,ZCCS,REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,Nedp,...
NONC,MXCCS,RMYU0,NE,CCS,ipconst)
RNOR = 0;
ZNOR = 0;
BR = zeros(1,NAPB);
BZ = zeros(1,NAPB);
PSIFLX = zeros(1,NFLX);
PSIC = zeros(1,MXCCS);
FC = zeros(1,n100);
%
%	
for I = 1:NAPB+NFLX	
	fprintf(WAHAHA,'%d Type=%d RS/ZS= %d %d Tet= %d\r\n',I,ITYPE(I),RS(I),ZS(I),TET(I));
end
%
% **********************************************************************
if (IT > 1)
else
% **********************************************************************
%����Ō��̈ʒu�ɖ߂�
for I = 1:CCS
    RCCS(I,NCCS(I)+1) = RCCS(I,1);
    ZCCS(I,NCCS(I)+1) = ZCCS(I,1);
end
%    !======================================================================
%    !
%    !    FFFFFF
%    !    FF                             | T-PROBE |
%    !    FFFF                VECTOR  FF=| N-PROBE |
%    !    FF                             |FLUX-LOOP|
%    !    FF                             |   CCS   |
%    !
%
    fprintf(WAHAHA,'%s\r\n','*****�Z���T�[�M������R�C���d����^�̍�����*****');
%  �����̉�FF���쐬�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I
%    !FLUX-LOOP 
    fprintf(WAHAHA,'%s\r\n','FLUX-LOOP     PSI  caused by external coils');
    for L = 1:NFLX
        PSIFLX(L) = 0.0D0;
        PSIFLX(L) = GETA;
        [PPSIFLX(1:KCMX),PHIR,PHIZ,PPSIA,PPSIB,PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,...
        DAG,DBG,DAH,DBH] = STARB(0,RS(L),ZS(L),RC(1:KCMX),ZC(1:KCMX),RNOR,ZNOR); % OK
        % �R�C����^�� = ���i��* * I_coil * ��0�j
	    PSIFLX(L) = PSIFLX(L) + sum(PPSIFLX(1:KCMX).*ECI(1:KCMX)*RMYU0);
        fprintf(WAHAHA,'%d %d\r\n',L,PSIFLX(L));
        FF_sub_fl(L+NAPB) = FF(L+NAPB);
        FF(L+NAPB) = FF(L+NAPB) - PSIFLX(L);%   ! ���ʏ������܂�
        FC(L+NAPB) = PSIFLX(L); %  ! �R�C���d����^
    end
    figure('Name','Flux function minus coil current contribution','NumberTitle','off')
    plot(FF_sub_fl(NAPB + 1 : NAPB + NFLX), 'DisplayName', 'Measured flux function')
    hold on
    plot(FF(NAPB + 1 : NAPB + NFLX), 'DisplayName', 'Measured flux function - calculated flux function by coil current')
    hold on
    plot(FC(NAPB + 1 : NAPB + NFLX), 'DisplayName', 'Calculated flux function by coil current')
    xlabel('Channel')
    ylabel('Flux function [Wb/rad]')
    legend
    
%    !T-PROBE & N-PROBE
    fprintf(WAHAHA,'%s\r\n','T-PROBE & N-PROBE   B  caused by external coils');
    for L = 1:NAPB
        BR(L) = 0.0D0;
        BZ(L) = 0.0D0;
        [PHI,PHIR,PHIZ,PHIA(1:KCMX),PHIB(1:KCMX),PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,DAG,DBG,...
        DAH,DBH] = STARB(1,RS(L+NFLX),ZS(L+NFLX),RC(1:KCMX),ZC(1:KCMX),RNOR,ZNOR); % OK
        % �R�C����^BR = ���i�i-�݃�*/��b / R�j * I_coil * ��0�j�@��b�F�Z���Tz�ʒu
        % �R�C����^BZ = ���i�i�݃�*/��a / R�j * I_coil * ��0�j�@��a�F�Z���Tr�ʒu
        BR(L) = sum((-PHIB(1:KCMX)/RS(L+NFLX)).*ECI(1:KCMX)*RMYU0);
        BZ(L) = sum((PHIA(1:KCMX)/RS(L+NFLX)).*ECI(1:KCMX)*RMYU0);
        % �v�����������̎�����Z�o
        BBB = BR(L)*cos(TET(L+NFLX))+BZ(L)*sin(TET(L+NFLX));
        fprintf(WAHAHA,'%d %d\r\n',L,BBB);
        FF_sub_f(L) = FF(L);
        FF(L) = FF(L) - BBB;             
        FC(L) = BBB; %   ! �R�C���d����^       
    end
    figure('Name','B minus coil current contribution','NumberTitle','off')
    plot(FF_sub_f(1 : NAPB), 'DisplayName', 'Measured axial field')
    hold on
    plot(FF(1 : NAPB), 'DisplayName', 'Measured axial field - calculated axial field by coil current')
    hold on
    plot(FC(1 : NAPB), 'DisplayName', 'Calculated axial field by coil current')
    xlabel('Channel')
    ylabel('Axial field [T]')
    legend
    
    
%    !CCS
    fprintf(WAHAHA,'%s\r\n','CCS       PSI  caused by external coils');
    for III = 1:CCS
        for L = 1:NCCN(III)
            PSIC(L) = GETA;
                [PPSIC(1:KCMX),PHIR,PHIZ,PPSIA,PPSIB,PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,DAG,...
            DBG,DAH,DBH] = STARB(0,RCCN(III,L),ZCCN(III,L),RC(1:KCMX),ZC(1:KCMX),RNOR,ZNOR); % OK
            % �R�C����^�� = ���i��* * I_coil * ��0�j
            PSIC(L) = PSIC(L) + sum(PPSIC(1:KCMX).*ECI(1:KCMX)*RMYU0);
            FF(L+sum(NCCN(1:III-1))+NAPB+NFLX) = - PSIC(L);%  !  ���ʏ������܂ށ@�@
            FC(L+sum(NCCN(1:III-1))+NAPB+NFLX) = PSIC(L);  %  ! �R�C���d����^�@�@
            fprintf(WAHAHA,'%d %d %d\r\n',L,PSIC(L),FF(L+sum(NCCN(1:III-1))+NAPB+NFLX));
        end
    end
%%
%!======================================================================
%!
%!    AA
%!   AAAA
%!  AA  AA                  |GT -HT|       | T-PROBE |
%! AA    AA       MATRIX AA=|GN -HN|  <--- | N-PROBE |
%! AAAAAAAA                 |GF -HF|       |FLUX-LOOP|
%! AA    AA                 |GC -HC|       |   CCS   |
%! AA    AA 
%!
%%%
    fid99 = fopen('MINDIST.txt','w'); %99
    frewind(fid99);
    fid100 = fopen('SEKIBUNCHECK.PRI','w'); %100 
    frewind(fid100);
    fprintf(fid99,'%s\r\n','****************************************************');
    fprintf(fid99,'%s\r\n','***    In the Subr. FORM ***************************');
    fprintf(fid99,'%s\r\n\r\n','****************************************************');
    fprintf(fid100,'%s\r\n','***************************************************');
    fprintf(fid100,'%s\r\n','***    In the Subr. FORM **************************');
    fprintf(fid100,'%s\r\n\r\n','***************************************************');
%%%%
%!FLUX-LOOP
%disp(NFLX)
% 
for III = 1:CCS
    for I = 1:NFLX
        for K = 1:NE(III)
            [HW,GW,GR,GZ,HR,HZ] = INTEGS(RS(I),ZS(I),RCCS(III,2*K-1),ZCCS(III,2*K-1),...
            RCCS(III,2*K),ZCCS(III,2*K),RCCS(III,2*K+1),ZCCS(III,2*K+1));% OK
            for JJ = 1:3
                KK = 3*(K-1)+JJ;
                %disp(GW(JJ))
                AA(I+NAPB,KK+3*sum(NE(1:III-1))) = AA(I+NAPB,KK+3*sum(NE(1:III-1)))+GW(JJ);
                AA(I+NAPB,KK+3*sum(NE(1:III-1))+sum(NCCN)) = AA(I+NAPB,KK+3*sum(NE(1:III-1))+sum(NCCN))-HW(JJ);
            end
        end% 119
    end
%C      
%!T-PROBE & N-PROBE
    for I = 1:NAPB
        % UTST�ł�TET=pi/2�iBz�݂̂��v�����Ă��邽�߁j
        COST=cos(TET(I+NFLX));
	    SINT=sin(TET(I+NFLX));
        for K = 1:NE(III)
            [HW,GW,GR,GZ,HR,HZ] = INTEGS(RS(I+NFLX),ZS(I+NFLX),RCCS(III,2*K-1),...
            ZCCS(III,2*K-1),RCCS(III,2*K),ZCCS(III,2*K),RCCS(III,2*K+1),ZCCS(III,2*K+1)); % OK
            for JJ = 1:3
                KK = 3*(K-1)+JJ;
                G = -COST*GZ(JJ)/RS(I+NFLX)+SINT*GR(JJ)/RS(I+NFLX);
                H = -COST*HZ(JJ)/RS(I+NFLX)+SINT*HR(JJ)/RS(I+NFLX);
                AA(I,KK+3*sum(NE(1:III-1))) = AA(I,KK+3*sum(NE(1:III-1)))+G;
                AA(I,KK+3*sum(NE(1:III-1))+sum(NCCN)) = AA(I,KK+3*sum(NE(1:III-1))+sum(NCCN))-H;
            end
        end%29 
    end   
%
    %!CCS
    for I = 1:NCCN(III)
        HII=0.0;
        for K = 1:NE(III)
            %���ٓ_���܂ރZ�O�����g�̋��E�ϕ��͓��ْl�ϕ������s
            if and((3*K) >= I , I >= (3*K-2));
                %//// ���ْl�ϕ� ///////////////////////////////////////////////////////
                if (I == (3*K-2)) % ���E��1�Ԗڂ̐ߓ_
                    NODO = 1;
                else
%                    if(I == (3*K-1))
	                    NODO = 2.*(I == (3*K-1)) +3.*(I ~= (3*K-1));% ���E��2�A3�Ԗڂ̐ߓ_
%                         disp(NODO)
%                         
%                    else
%	                    NODO = 3;
%                    end
                end
                [GW,HW] = INLOGSA(RCCN(III,I),ZCCN(III,I),RCCS(III,2*K-1),ZCCS(III,2*K-1),RCCS(III,2*K),...
                ZCCS(III,2*K),RCCS(III,2*K+1),ZCCS(III,2*K+1),NODO); % OK
            else
                %//// �ʏ�̐ϕ� ///////////////////////////////////////////////////////
                [HW,GW,GR,GZ,HR,HZ] = INTEGS(RCCN(III,I),ZCCN(III,I),RCCS(III,2*K-1),ZCCS(III,2*K-1),...
                RCCS(III,2*K),ZCCS(III,2*K),RCCS(III,2*K+1),ZCCS(III,2*K+1)); % OK
            end
            for JJ = 1:3
	            KK = 3*(K-1)+JJ;
                AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,KK+3*sum(NE(1:III-1))) = AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,KK+3*sum(NE(1:III-1)))+GW(JJ);
                AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,KK+3*sum(NE(1:III-1))+sum(NCCN)) = AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,KK+3*sum(NE(1:III-1))+sum(NCCN))-HW(JJ);
            end
        end
        for J = 1:NCCN(III)
            if (J == I)
                continue
            end
	        HII=HII+AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,J+sum(NCCN(1:III-1))+sum(NCCN));
        end% 101
%        if (HII < -0.001)
            HII = HII + 1.0.*(HII < -0.001);
%        end
        AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,I+sum(NCCN(1:III-1))+sum(NCCN)) = -HII;
    end % 140
end

    % ??????????????????????????????????????????????????????????????????    
    %  �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��
    %    �@�^��e��@�@�@�@�@�^��e��@�@�@�@�@�^��e��@�@�@�@�^��e��
    %     KNE=No. of boundary elements along the vauum vessel
    %     KNN=No. of nodes along the vauum vessel (KNN=KNE*2)
    %     (REV(),ZEV())=Eddy Current Nodes on the vacuum vessel
    % ??????????????????????????????????????????????????????????????????    
    %****
    JJJJ = sum(NCCN)+sum(NCCN);
    if (KNE <= 0) %! GOTO 991
    else
        AMYU0 = RMYU0*1.0D06;   %! NAMUAMUdabutsu  #1
        %****
        %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
        %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
        %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
        %    ��K��(Non Conforming)�Q�d���v�f�̎���  (if NONC=1)
        if (NONC == 0) %GOTO 990
            fprintf('%s\r\n','�^��e���̉Q�d���� flux loop �Z���T�[�ɍ�郵');
            for I = 1:NFLX
                A = RS(I);
                B=ZS(I);
                for K = 1:2:KNN-1
                [GW,GR,GZ] =  EXTINDC(A,B,REV(K),ZEV(K),REV(K+1),ZEV(K+1),REV(K+2),...
                ZEV(K+2),NONC,fid99,fid100); % OK
                    for JJ = 1:3
                        EE = GW(JJ)*AMYU0;     % ! flux*AMYU0         
                        if and(K == KNN-1, JJ == 3)
                            AA(I+NAPB,sum(NCCN)*2+1) = AA(I+NAPB,sum(NCCN)*2+1)+EE;    %        ! ### 2
                            AA(I+NAPB,sum(NCCN)*2+K-1+JJ) = AA(I+NAPB,sum(NCCN)*2+K-1+JJ)+EE;%  ! ### 2
                        else
                            
                        end
                    end
                end
            end
%
            fprintf('%s\r\n','�^��e���̉Q�d���� ����Z���T�[�ɍ��a');
            for I = 1:NAPB
                COST = cos(TET(I+NFLX));
                SINT = sin(TET(I+NFLX));
                A = RS(I+NFLX);
                B = ZS(I+NFLX);
                for K = 1:2:KNN-1
                    [GW,GR,GZ] =  EXTINDC(A,B,REV(K),ZEV(K),REV(K+1),ZEV(K+1),...
                    REV(K+2),ZEV(K+2),NONC,fid99,fid100); % OK
                    for JJ = 1:3
                        EE = (-COST*GZ(JJ)+SINT*GR(JJ))*AMYU0/A;   %! AMYU0/A
%%!!   (-COST,SINT) ------ Need to reconfirm --- OK!!
                        if and(K == KNN-1, JJ == 3)
                            AA(I,sum(NCCN)*2+1) = AA(I,sum(NCCN)*2+1)+EE;     %   ! ### 1
                        else
                            AA(I,sum(NCCN)*2+K-1+JJ) = AA(I,sum(NCCN)*2+K-1+JJ)+EE; % ! ### 1
                        end
                    end
                end
            end
%
            fprintf('%s\r\n','�^��e���̉Q�d���� CCS�ɍ�郵');
            for III = 1:CCS
                for I = 1:NCCN(III)
                    A = RCCN(III,I);
                    B = ZCCN(III,I);
                    for K = 1:2:KNN-1
                        [GW,GR,GZ] =  EXTINDC(A,B,REV(K),ZEV(K),REV(K+1),ZEV(K+1),...
                        REV(K+2),ZEV(K+2),NONC,fid99,fid100); % OK
                        for JJ = 1:3
                            EE = GW(JJ)*AMYU0;    % ! flux*AMYU0
                            if and(K == KNN-1, JJ == 3)
                                AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+1) = AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+1)+EE;            %! ### 3
                            else
                                AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+K-1+JJ) = AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+K-1+JJ)+EE; % ! ### 3
                            end
                        end
                    end
                end 
            end
            JJJJ = sum(NCCN)+sum(NCCN)+KNN;
%
        else
            fprintf('%s\r\n','�^��e���̉Q�d���� flux loop �Z���T�[�ɍ�郵');
            for I = 1:NFLX
                A=RS(I);
                B=ZS(I);
                for K = 1:KNE
                    [GW,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),REV(2*K+1),...
                    ZEV(2*K+1),NONC,fid99,fid100); % OK
%                     if (abs(REV(2*K) -0.695) < 0.001 )
%                         [GW_1,GR,GZ] =  EXTINDC(A,B,REV(2*K-1)-0.002,ZEV(2*K-1),REV(2*K)-0.002,ZEV(2*K),REV(2*K+1)-0.002,...
%                         ZEV(2*K+1),NONC,fid99,fid100); % OK                         
%                         [GW_2,GR,GZ] =  EXTINDC(A,B,REV(2*K-1)+0.002,ZEV(2*K-1),REV(2*K)+0.002,ZEV(2*K),REV(2*K+1)+0.002,...
%                         ZEV(2*K+1),NONC,fid99,fid100); % OK
%                         GW = GW + GW_1 + GW_2;
%                     end
%                     if (abs(abs(ZEV(2*K)) -0.2925) < 0.005 )
%                         [GW_1,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.00625,REV(2*K),ZEV(2*K)+0.00625,REV(2*K+1),...
%                         ZEV(2*K+1)+0.00625,NONC,fid99,fid100); % OK                         
%                         [GW_2,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)-0.00625,REV(2*K),ZEV(2*K)-0.00625,REV(2*K+1),...
%                         ZEV(2*K+1)-0.00625,NONC,fid99,fid100); % OK
%                         GW = GW + GW_1 + GW_2;
%                     end
%                     if (abs(REV(2*K) -0.1183365) < 0.001)
%                         if(ZEV(2*K) > 0)
%                             fugou = 1;
%                         else
%                             fugou = -1;
%                         end
%                             [GW_1,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.042*fugou,REV(2*K),ZEV(2*K)+0.042*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.042*fugou,NONC,fid99,fid100); % OK                         
%                             [GW_2,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.084*fugou,REV(2*K),ZEV(2*K)+0.084*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.084*fugou,NONC,fid99,fid100); % OK
%                             [GW_3,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.126*fugou,REV(2*K),ZEV(2*K)+0.126*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.126*fugou,NONC,fid99,fid100); % OK                         
%                             [GW_4,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.168*fugou,REV(2*K),ZEV(2*K)+0.168*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.168*fugou,NONC,fid99,fid100); % OK
%                         GW = GW + GW_1 + GW_2 + GW_3 + GW_4;
%                     end
%                     if (abs(ZEV(2*K)) < 0.001 )
%                         GW = GW*0;
%                     end
                    for JJ=1:3
                        KK=3*(K-1)+JJ;
                        EE = GW(JJ)*AMYU0;%      ! flux*AMYU0         
                        AA(I+NAPB,sum(NCCN)*2+KK)=AA(I+NAPB,sum(NCCN)*2+KK)+EE;%  ! ### 2
%                         fprintf('kondo')
%                         I+NAPB
%                         sum(NCCN)*2+KK
%                         
                    end
                end
            end
%
            fprintf('%s\r\n','�^��e���̉Q�d���� ����Z���T�[�ɍ��a');
            for I = 1:NAPB
                COST = cos(TET(I+NFLX));
                SINT = sin(TET(I+NFLX));
                A = RS(I+NFLX);
                B = ZS(I+NFLX);
                for K = 1:KNE
                    [GW,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),REV(2*K+1),...
                    ZEV(2*K+1),NONC,fid99,fid100); % OK
%                     if (abs(REV(2*K) -0.695) < 0.001 )
%                         [GW,GR_1,GZ_1] =  EXTINDC(A,B,REV(2*K-1)-0.002,ZEV(2*K-1),REV(2*K)-0.002,ZEV(2*K),REV(2*K+1)-0.002,...
%                         ZEV(2*K+1),NONC,fid99,fid100); % OK                         
%                         [GW,GR_2,GZ_2] =  EXTINDC(A,B,REV(2*K-1)+0.002,ZEV(2*K-1),REV(2*K)+0.002,ZEV(2*K),REV(2*K+1)+0.002,...
%                         ZEV(2*K+1),NONC,fid99,fid100); % OK
%                         GR = GR + GR_1 + GR_2;
%                         GZ = GZ + GZ_1 + GZ_2;
%                     end
%                     if (abs(abs(ZEV(2*K)) -0.2925) < 0.005 )
%                         [GW,GR_1,GZ_1] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.00625,REV(2*K),ZEV(2*K)+0.00625,REV(2*K+1),...
%                         ZEV(2*K+1)+0.00625,NONC,fid99,fid100); % OK                         
%                         [GW,GR_2,GZ_2] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)-0.00625,REV(2*K),ZEV(2*K)-0.00625,REV(2*K+1),...
%                         ZEV(2*K+1)-0.00625,NONC,fid99,fid100); % OK
%                         GR = GR + GR_1 + GR_2;
%                         GZ = GZ + GZ_1 + GZ_2;
%                     end
%                     if (abs(REV(2*K) -0.1183365) < 0.001)
%                         if(ZEV(2*K) > 0)
%                             fugou = 1;
%                         else
%                             fugou = -1;
%                         end
%                             [GW_1,GR_1,GZ_1] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.042*fugou,REV(2*K),ZEV(2*K)+0.042*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.042*fugou,NONC,fid99,fid100); % OK                         
%                             [GW_2,GR_2,GZ_2] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.084*fugou,REV(2*K),ZEV(2*K)+0.084*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.084*fugou,NONC,fid99,fid100); % OK
%                             [GW_3,GR_3,GZ_3] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.126*fugou,REV(2*K),ZEV(2*K)+0.126*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.126*fugou,NONC,fid99,fid100); % OK                         
%                             [GW_4,GR_4,GZ_4] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.168*fugou,REV(2*K),ZEV(2*K)+0.168*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.168*fugou,NONC,fid99,fid100); % OK
%                         GR = GR + GR_1 + GR_2 + GR_3 + GR_4;
%                         GZ = GZ + GZ_1 + GZ_2 + GZ_3 + GZ_4;
%                     end
%                     if (abs(ZEV(2*K)) < 0.001 )
%                         GR = GR*0;
%                         GZ = GZ*0;
%                     end
                    for JJ = 1:3
                        KK=3*(K-1)+JJ;
                        EE = (-COST*GZ(JJ)+SINT*GR(JJ))*AMYU0/A; %  ! AMYU0/A
                        %!!   (-COST,SINT) ------ Need to reconfirm --- OK!!
                        AA(I,sum(NCCN)*2+KK)=AA(I,sum(NCCN)*2+KK)+EE;%  ! ### 1
                    end
                end
            end
%
            fprintf('%s\r\n','�^��e���̉Q�d���� CCS�ɍ�郵');
            for III = 1:CCS
                for I = 1:NCCN(III)
                    A = RCCN(III,I);
                    B = ZCCN(III,I);
                    for K = 1:KNE
                        [GW,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),REV(2*K+1),...
                        ZEV(2*K+1),NONC,fid99,fid100); % OK
%                         if (abs(REV(2*K) -0.695) < 0.001 )
%                             [GW_1,GR,GZ] =  EXTINDC(A,B,REV(2*K-1)-0.002,ZEV(2*K-1),REV(2*K)-0.002,ZEV(2*K),REV(2*K+1)-0.002,...
%                             ZEV(2*K+1),NONC,fid99,fid100); % OK                         
%                             [GW_2,GR,GZ] =  EXTINDC(A,B,REV(2*K-1)+0.002,ZEV(2*K-1),REV(2*K)+0.002,ZEV(2*K),REV(2*K+1)+0.002,...
%                             ZEV(2*K+1),NONC,fid99,fid100); % OK
%                             GW = GW + GW_1 + GW_2;
%                         end
%                         if (abs(abs(ZEV(2*K)) -0.2925) < 0.005 )
%                             [GW_1,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.00625,REV(2*K),ZEV(2*K)+0.00625,REV(2*K+1),...
%                             ZEV(2*K+1)+0.00625,NONC,fid99,fid100); % OK                         
%                             [GW_2,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)-0.00625,REV(2*K),ZEV(2*K)-0.00625,REV(2*K+1),...
%                             ZEV(2*K+1)-0.00625,NONC,fid99,fid100); % OK
%                             GW = GW + GW_1 + GW_2;
%                         end
%                         if (abs(REV(2*K) -0.1183365) < 0.001)
%                             if(ZEV(2*K) > 0)
%                                 fugou = 1;
%                             else
%                                 fugou = -1;
%                             end
%                             [GW_1,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.042*fugou,REV(2*K),ZEV(2*K)+0.042*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.042*fugou,NONC,fid99,fid100); % OK                         
%                             [GW_2,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.084*fugou,REV(2*K),ZEV(2*K)+0.084*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.084*fugou,NONC,fid99,fid100); % OK
%                             [GW_3,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.126*fugou,REV(2*K),ZEV(2*K)+0.126*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.126*fugou,NONC,fid99,fid100); % OK                         
%                             [GW_4,GR,GZ] =  EXTINDC(A,B,REV(2*K-1),ZEV(2*K-1)+0.168*fugou,REV(2*K),ZEV(2*K)+0.168*fugou,REV(2*K+1),...
%                             ZEV(2*K+1)+0.168*fugou,NONC,fid99,fid100); % OK
%                             GW = GW + GW_1 + GW_2 + GW_3 + GW_4;
%                        end
%                         if (abs(ZEV(2*K)) < 0.001 )
%                             GW = GW*0;
%                         end
                        for JJ = 1:3
                            KK = 3*(K-1)+JJ;
                            EE = GW(JJ)*AMYU0;%     ! flux*AMYU0
                            AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+KK)=AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+KK)+EE;%  ! ### 3
                        end
                    end
                end
            end
            JJJJ=sum(NCCN)+sum(NCCN)+KNE*3;%             
        end
    end
 
    
% disp(size(FF))
% figure('Name','Flux Profile','NumberTitle','off')
% plot(FF)

    
disp(size(AA))
% figure('Name','Flux Profile','NumberTitle','off')
% contourf(AA,100); 
% c = colorbar;
% caxis([-0.01 0.01]);


disp(AA(78,:))

disp(AA(79,:))

    

%�@�֌W�Ȃ��悤�Ȃ̂ŃR�����g�A�E�g�i�ߓ��j
    % ??????????????????????????????????????????????????????????????????    
    %  �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��
    %    �@���艻�@�@�@�@�@���艻�@�@�@�@�@���艻�@�@�@�@���艻��
    %CAUTION!! The stabilizer is not closed in the poloidal direction.
    %     KSE=No. of boundary elements along the stabilizer
    %     KSN=No. of nodes along the stabilizer (KSN=KNE*2+1)
    %     (RES(),ZES())=Eddy Current Nodes on the stabilizer
    % ??????????????????????????????????????????????????????????????????    
%     if (KSE <= 0) % 991 to 992
%     else
%         AMYU0 = RMYU0*1.0D06;  % ! NAMUAMUdabutsu  #2
%         fprintf('%s\r\n','���艻��̉Q�d���� flux loop �Z���T�[�ɍ�郵');
%         for I = 1:NFLX
%             A = RS(I);
%             B = ZS(I);
%             for K = 1:2:KSN-2
%                 [GW,GR,GZ] =  EXTINDC(A,B,RES(K),ZES(K),RES(K+1),ZES(K+1),RES(K+2),...
%                 ZES(K+2),NONC,fid99,fid100);% OK 
%                 for JJ = 1:3
%                     EE = GW(JJ)*AMYU0;    %  ! flux*AMYU0         
%                     AA(I+NAPB,sum(NCCN)*2+KNN+K-1+JJ) = AA(I+NAPB,sum(NCCN)*2+KNN+K-1+JJ)+EE; % ! ### 2
%                 end
%             end
%         end
% %CC
%         fprintf('%s\r\n','���艻��̉Q�d���� ����Z���T�[�ɍ��a');
%         for I = 1:NAPB
%             COST = cos(TET(I+NFLX));
%             SINT = sin(TET(I+NFLX));
%             A = RS(I+NFLX);
%             B = ZS(I+NFLX);
%             for K = 1:2:KSN-2
%                 [GW,GR,GZ] =  EXTINDC(A,B,RES(K),ZES(K),RES(K+1),ZES(K+1),RES(K+2),...
%                 ZES(K+2),NONC,fid99,fid100); % OK 
%                 for JJ = 1:3
%                     EE = (-COST*GZ(JJ)+SINT*GR(JJ))*AMYU0/A;  % ! AMYU0/A
%                     AA(I,sum(NCCN)*2+KNN+K-1+JJ) = AA(I,sum(NCCN)*2+KNN+K-1+JJ)+EE;  %! ### 1
%                 end
%             end
%         end
% %CC
%         fprintf('%s\r\n','���艻��̉Q�d���� CCS�ɍ�郵');
%         for III = 1:CCS
%            for I = 1:NCCN(III)
%                A = RCCN(III,I);
%                B = ZCCN(III,I);
%                for K = 1:2:KSN-2
%                    [GW,GR,GZ] =  EXTINDC(A,B,RES(K),ZES(K),RES(K+1),ZES(K+1),RES(K+2),...
%                    ZES(K+2),NONC,fid99,fid100); % OK
%                    for JJ = 1:3
%                        EE = GW(JJ)*AMYU0;   %  ! flux*AMYU0
%                        AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+KNN+K-1+JJ) = AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+KNN+K-1+JJ)+EE; % ! ### 3
%                    end
%                end 
%            end  
%        end
%        JJJJ = sum(NCCN)+sum(NCCN)+sum(KNN)+sum(KSN);
% %
% % ??????????????????????????????????????????????????????????????????    
% %
%          
%     end
% TOTAL IP (USHIKI)
if (ipconst == 1)
%FF(sum(NCCN)+NAPB+NFLX+1) = 50000*RMYU0;
    for III = 1:CCS
        for K = 1:NE(III)
            [HIP] = INTEIP(RCCS(III,2*K-1),ZCCS(III,2*K-1),...
            RCCS(III,2*K),ZCCS(III,2*K),RCCS(III,2*K+1),ZCCS(III,2*K+1));% OK
            for JJ = 1:3
                KK = 3*(K-1)+JJ;
                AA(sum(NCCN)+NAPB+NFLX+1,KK+3*sum(NE(1:III-1))) = AA(sum(NCCN)+NAPB+NFLX+1,KK+3*sum(NE(1:III-1))) + HIP(JJ);
                AA(sum(NCCN)+NAPB+NFLX+1,KK+3*sum(NE(1:III-1))+sum(NCCN)) = 0;
            end
        end% 119
    end
end
% TOTAL IP(USHIKI)
%C      WRITE(IPR,16) NFLX+NAPB+NCCS
if (ipconst == 1)
    fprintf(WAHAHA,'No. of information data = %d\r\n',NFLX+NAPB+sum(NCCN)+1);% +1�͑��d��
    fprintf(WAHAHA,'Vector FF(I) at the initial stage of iteration\r\n');
    for I = 1:NFLX+NAPB+sum(NCCN)+1
        fprintf(WAHAHA,'%d\r\n',FF(I));
    end
    fprintf(WAHAHA,'No. of unknowns =   %d   %d   %d   %d\r\n\r\n',JJJJ,sum(NCCN),sum(KNN),sum(KSN));
    fprintf(WAHAHA,'Matrix AA(I,J)\r\n');
    for I = 1:NFLX+NAPB+sum(NCCN)+1
 	    fprintf(WAHAHA,'I=%d\r\n',I);
        for J = JJJJ
            fprintf(WAHAHA,'%d\r\n',AA(I,J));
        end
    end
else    
%C      WRITE(IPR,16) NFLX+NAPB+NCCS
    fprintf(WAHAHA,'No. of information data = %d\r\n',NFLX+NAPB+sum(NCCN));
    fprintf(WAHAHA,'Vector FF(I) at the initial stage of iteration\r\n');
    for I = 1:NFLX+NAPB+sum(NCCN)
        fprintf(WAHAHA,'%d\r\n',FF(I));
    end
    fprintf(WAHAHA,'No. of unknowns =   %d   %d   %d   %d\r\n\r\n',JJJJ,sum(NCCN),KNN,KSN);
    fprintf(WAHAHA,'Matrix AA(I,J)\r\n');
    for I = 1:NFLX+NAPB+sum(NCCN)
 	    fprintf(WAHAHA,'I=%d\r\n',I);
        for J = JJJJ
            fprintf(WAHAHA,'%d\r\n',AA(I,J));
        end
    end
end
% ??????????????????????????????????????????????????????????????????    
%
end
FF(1+NAPB:NFLX+NAPB) = FF(1+NAPB:NFLX+NAPB)-GETA+OLDGT;    %! �v�Ċm�F
FF(1+NAPB+NFLX:sum(NCCN)+NAPB+NFLX) = FF(1+NAPB+NFLX:sum(NCCN)+NAPB+NFLX)-GETA+OLDGT;%    ! �v�Ċm�F
end
%
%% INTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [PSI,DELGE,RCCS,ZCCS,XPSI] = INTER(IGOAL,GETA,RCCS,ZCCS,CR,CZ,FFOUT,...
n100,n50,RC,ZC,ECIGRP,ECI,KCMX,NAPB,NFLX,FLXLP,BSNSR,RS,ZS,TET,RCCN,ZCCN,NCCN,...
KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,WAHAHA,MXCCS,MXINT,...
NCCS,NINT,RMYU0,NE,CCS)
%%
PSI = zeros(1,MXINT);
PSIA = zeros(1,MXINT);
PSIB = zeros(1,MXINT);
RCCSR = zeros(1,MXCCS+1);
ZCCSR = zeros(1,MXCCS+1);
FI = zeros(1,MXCCS);
DFI = zeros(1,MXCCS);%  ! BOUNDARY CONDITION FI:��, DFI:d��/dn
XPSI = zeros(1,n100);
XBBR = zeros(1,2000);
XBBZ = zeros(1,2000);
DELGE = 0; % ushiki
%
    fprintf(fid99,'%s\r\n','/'); 
    fprintf(fid100,'%s\r\n','/'); 
    fprintf(fid99,'%s\r\n', '****************************************************');
    fprintf(fid99,'%s\r\n', '***    In the Subr. INTER **************************');
    fprintf(fid99,'%s\r\n', '****************************************************');
    fprintf(fid99,'%s\r\n','/'); 
    fprintf(fid100,'%s\r\n', '****************************************************');
    fprintf(fid100,'%s\r\n', '***    In the Subr. INTER **************************');
    fprintf(fid100,'%s\r\n', '****************************************************');
    fprintf(fid100,'%s\r\n','/'); 
%
    for I = 1:CCS
        RCCS(I,NCCS(I)+1) = RCCS(I,1);
        ZCCS(I,NCCS(I)+1) = ZCCS(I,1);
    end
%
    DFI(1:sum(NCCN)) = FFOUT(1:sum(NCCN));
    FI(1:sum(NCCN)) = FFOUT(1 + sum(NCCN):sum(NCCN) + sum(NCCN));
    fprintf('%s %d\r\n', 'GETA in INTER =',GETA);
        if (IGOAL > 0) %GOTO 999
        else
           %%  ����Z���T�[�ɍ��a'')')
             for L = 1:NAPB
                 A = RS(NFLX+L);
                 B = ZS(NFLX+L);
                 % *******************************************************************
                 [PSIL,PSIA,PSIB] = QINTER(A,B,GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,...
                 n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,...
                 ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS); % OK
                 % ********************************************************************
                 XBBR(L) = -PSIB/A;
                 XBBZ(L) =  PSIA/B;
             end
%                  PSIA = zeros(NAPB); 
%                  PSIB = zeros(NAPB);
%                  L = 1:NAPB
% %                 A = RS(NFLX+L);
% %                 B = ZS(NFLX+L);
%                  % *******************************************************************
%                  [PSIL,PSIA(L),PSIB(L)] = QINTER(RS(NFLX+L),ZS(NFLX+L),GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,...
%                  n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,...
%                  ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE); % OK
%                  % ********************************************************************
%                  XBBR(L) = -PSIB(L)./RS(NFLX+L);
%                  XBBZ(L) =  PSIA(L)./ZS(NFLX+L);           
%     
	         fprintf(WAHAHA,'%s\r\n','Reproducibility of field sensor signals');
% 
	         fid104 = fopen('Comparison_TotalFieldSignal.txt','w');%104
	         fid106 = fopen('Discrepant_Field_Points.txt','w');%106
             II = 0;
             for I = 1:NAPB
 	             XBBB = XBBR(I)*cos(TET(I+NFLX))+XBBZ(I)*sin(TET(I+NFLX));
                 fprintf(WAHAHA, '%d Measured = %d Calculated = %d\r\n', I,BSNSR(I),XBBB);
                 fprintf(fid104,'%d %d\r\n',BSNSR(I),XBBB);
                 %�G���[
                 EE = 100.0D0*abs((XBBB-BSNSR(I))/BSNSR(I));
                 %�G���[��30���ȏ�̎�
                 if (EE > 30.0D0)
                     II = II+1;
                     fprintf(fid106,'%d %d\r\n', RS(I+NFLX),ZS(I+NFLX));
                     fprintf(WAHAHA,'%d %d %d %d %e\r\n', RS(I+NFLX),ZS(I+NFLX),BSNSR(I),XBBB,EE);
                 end
             end
             if (II <= 0)
                 SSS=1.0D10;
                 fprintf(fid106,'%d %d\r\n', SSS,SSS);
             end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
                %  Flux loop �Z���T�[�ɍ�郵'')')
             DELGE = 0.0D0;
             KNT = 0;
             for L = 1:NFLX
                 KNT=KNT+1;
                 A=RS(L);
                 B=ZS(L);
% *******************************************************************
                 [XPSI(L),PSIA,PSIB] = QINTER(A,B,GETA,RCCS,ZCCS,FFOUT,FI,...
                 DFI,n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,...
                 KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS); % OK
% ********************************************************************CC 
                 DELGE = DELGE+FLXLP(L)-XPSI(L);
             end
%              L = 1:NFLX;
%              KNT = L;
% %             A=RS(L);
% %                 B=ZS(L);
% % *******************************************************************
%                  [XPSI(L),PSIA,PSIB] = QINTER(RS(L),ZS(L),GETA,RCCS,ZCCS,FFOUT,FI,...
%                  DFI,n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,...
%                  KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE); % OK
% % ********************************************************************CC 
%                  DELGE = DELGE+sum(FLXLP(L)-XPSI(L));
%        
%
             CNT0 = KNT;
             DELGE = DELGE/CNT0;
             fid107 = fopen('Discrepant_Flux_Points.txt','w');%107
           	 fprintf(WAHAHA,'%s\r\n','Reproducibility of flux loop signals');
	         fid105 = fopen('Comparison_TotalFluxSignal.txt','w');%105
	         II = 0;
             for I = 1:NFLX
                 fprintf(WAHAHA, '%d Measured = %d Calculated = %d\r\n', I,FLXLP(I),XPSI(I));
                 fprintf(fid105,'%d %d\r\n', FLXLP(I),XPSI(I));
                 EE = 100.0D0*abs((XPSI(I)-FLXLP(I))/FLXLP(I));
                 if (EE > 30.0D0)
                     II = II+1;
                     fprintf(fid107,'%d %d\r\n',RS(I),ZS(I));
                     fprintf(WAHAHA,'%d %d %d %d %d\r\n', RS(I),ZS(I),FLXLP(I),XPSI(I),EE);
                 end
             end
             if (II <= 0)
                 SSS = 1.0D10;
                 fprintf(fid107,'%d %d\r\n', SSS,SSS);
             end
             if (IGOAL <= 0)
                 return
             end
         end
         % *****************************************************************
         % *****************************************************************
         %   �������z�}�b�v�̍쐬
         % *****************************************************************
         % *****************************************************************
         %  Sum of contributions by CCS & external coils
         %!INTER INTER INTER
         %  �R�C���d�������_�ɍ�郵
         for L = 1:NINT
             A = CR(L);
             B = CZ(L);
             % *******************************************************************
             [PSI(L),PSIA,PSIB] = QINTER(A,B,GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
             n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,...
             KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS); % OK
             % ******************************************************************* 
         end         
%              L = 1:NINT;
% %             A = CR(L);
% %             B = CZ(L);
%              % *******************************************************************
%              [PSI(L),PSIA,PSIB] = QINTER(CR(L),CZ(L),GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
%              n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,...
%              KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE); % OK
%              % ******************************************************************* 
% 
%          fid111 = fopen('CCSR_Check.txt','w');
% %          for I = 1:NCCS %   !RCCSR,ZCCSR�i�����v���j��RCCS,ZCCS(���v���)
% %              IR = NCCS+1-I;
% %              RCCSR(I) = RCCS(IR);
% %              ZCCSR(I) = ZCCS(IR);
% %              fprintf(fid111,'%d %d\r\n', RCCSR(I),ZCCSR(I));
% %          end
% %image(unique(CR),unique(CZ),PSI);
%      for III = 1:numel(NCCS)
%          I = 1:NCCS(III); %   !RCCSR,ZCCSR�i�����v���j��RCCS,ZCCS(���v���)
%          IR = NCCS(III)+1-I;
%          RCCSR(III,I) = RCCS(III,IR);
%          ZCCSR(III,I) = ZCCS(III,IR);
%          fprintf(fid111,'%d %d\r\n', horzcat(RCCSR(III,I)',ZCCSR(III,I)'));
%          fprintf(fid111,'%d %d\r\n', RCCSR(III,1),ZCCSR(III,1));
%          fprintf(fid111,'%s\r\n', '==');
%          fprintf(fid111,'%d %d %d\r\n', NINT,MXCCS,NCCS(III));
%          for I = 1:NINT
%              [WWW] = OUTIN(RCCSR,ZCCSR,CR(I),CZ(I),MXCCS,NCCS); % OK
%              if (WWW > 0.6D0)
%                  PSI(I) = FI(1);
%                  fprintf(fid111,'%s %d\r\n', 'WWW =',WWW);
%              end
%          end
%      end
end
%
%% QINTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [PSIS,PSISA,PSISB] = QINTER(AS,BS,GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,...
         n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,...
         Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS)
%   PSIS = GETA;
   PPSI  = zeros(1,KCMX);
   PPSIA = zeros(1,KCMX);
   PPSIB = zeros(1,KCMX);
%   PSISA = 0.0D0;
%   PSISB = 0.0D0;
   RNOR = 0; % ushiki
   ZNOR = 0; % ushiki
%    
   %% �R�C���d������鎥��/����̉��Z
   [PPSI(1:KCMX),PHIR,PHIZ,PPSIA(1:KCMX),PPSIB(1:KCMX),PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,DAG,DBG,DAH,DBH]...
   = STARB(1,AS,BS,RC(1:KCMX),ZC(1:KCMX),RNOR,ZNOR); % OK
   % �܂��̓R�C���d����^�����Z
   PSIS  = GETA + sum(PPSI(1:KCMX).*ECI(1:KCMX).*RMYU0);
%    PSIS  = sum(PPSI(1:KCMX).*ECI(1:KCMX).*RMYU0);
   PSISA = sum(PPSIA(1:KCMX).*ECI(1:KCMX).*RMYU0);
   PSISB =  sum(PPSIB(1:KCMX).*ECI(1:KCMX).*RMYU0);
   % ??????????????????????????????????????????????????????????????????    
   %% �Q�d����^�̉��Z (1) ��
   % ??????????????????????????????????????????????????????????????????    
   %    �^��e���̉Q�d�������
   if (KNE > 0)  
   %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   %    ��K��(Non Conforming)�Q�d���v�f�̎���  (if NONC=1)  #1
       if (NONC == 0)% GOTO 990
           for K = 1:KNE      
               [GW,GR,GZ] = EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),...
               REV(2*K+1),ZEV(2*K+1),NONC,fid99,fid100); % OK
               for JJ = 1:3
	               KJ2 = 2*K-2+JJ;
                   if (KJ2 > KNN)
                       KJ2 = JJ-2;
                   end
	               PSIS = PSIS + GW(JJ)*FFOUT(sum(NCCN)*2+KJ2)*AMYU0;
	               PSISA = PSISA+ GR(JJ)*FFOUT(sum(NCCN)*2+KJ2)*AMYU0;
	               PSISB = PSISB+ GZ(JJ)*FFOUT(sum(NCCN)*2+KJ2)*AMYU0;
               end
           end
%                K = 1:KNE;
%                [GW(1:3,K),GR(1:3,K),GZ(1:3,K)] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),...
%                REV(2*K+1),ZEV(2*K+1),NONC,fid99,fid100); % OK
%                JJ = 1:3;
% 	           KJ2(K,1) = 2.*K-2+1;
% 	           KJ2(K,2) = 2.*K-2+2;
% 	           KJ2(K,3) = 2.*K-2+3;
%                KJ2(K,1) = (1-2).*(KJ2(K,1) > KNN) + KJ2(K,1).*(KJ2(K,1) <= KNN);
%                KJ2(K,2) = (2-2).*(KJ2(K,1) > KNN) + KJ2(K,1).*(KJ2(K,1) <= KNN);
%                KJ2(K,3) = (3-2).*(KJ2(K,1) > KNN) + KJ2(K,1).*(KJ2(K,1) <= KNN);
% 	           PSIS = PSIS + sum(sum(GW(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
% 	           PSISA = PSISA + sum(sum(GR(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
% 	           PSISB = PSISB + sum(sum(GZ(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
         
       else
           for K = 1:KNE
               [GW,GR,GZ] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),...
               REV(2*K+1),ZEV(2*K+1),NONC,fid99,fid100); % OK
%                if (abs(REV(2*K) -0.695) < 0.001 )
%                    [GW_1,GR_1,GZ_1] =  EXTINDC(AS,BS,REV(2*K-1)-0.002,ZEV(2*K-1),REV(2*K)-0.002,ZEV(2*K),REV(2*K+1)-0.002,...
%                    ZEV(2*K+1),NONC,fid99,fid100); % OK                         
%                    [GW_2,GR_2,GZ_2] =  EXTINDC(AS,BS,REV(2*K-1)+0.002,ZEV(2*K-1),REV(2*K)+0.002,ZEV(2*K),REV(2*K+1)+0.002,...
%                    ZEV(2*K+1),NONC,fid99,fid100); % OK
%                    GW = GW + GW_1 + GW_2;
%                    GR = GR + GR_1 + GR_2;
%                    GZ = GZ + GZ_1 + GZ_2;
%                end
%                if (abs(abs(ZEV(2*K)) -0.2925) < 0.005 )
%                    [GW_1,GR_1,GZ_1] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.00625,REV(2*K),ZEV(2*K)+0.00625,REV(2*K+1),...
%                    ZEV(2*K+1)+0.00625,NONC,fid99,fid100); % OK                         
%                    [GW_2,GR_2,GZ_2] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)-0.00625,REV(2*K),ZEV(2*K)-0.00625,REV(2*K+1),...
%                    ZEV(2*K+1)-0.00625,NONC,fid99,fid100); % OK
%                    GW = GW + GW_1 + GW_2;
%                    GR = GR + GR_1 + GR_2;
%                    GZ = GZ + GZ_1 + GZ_2;
%                end
%                if (abs(REV(2*K) -0.1183365) < 0.001)
%                    if(ZEV(2*K) > 0)
%                        fugou = 1;
%                    else
%                        fugou = -1;
%                    end
%                    [GW_1,GR_1,GZ_1] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.042*fugou,REV(2*K),ZEV(2*K)+0.042*fugou,REV(2*K+1),...
%                    ZEV(2*K+1)+0.042*fugou,NONC,fid99,fid100); % OK                         
%                    [GW_2,GR_2,GZ_2] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.084*fugou,REV(2*K),ZEV(2*K)+0.084*fugou,REV(2*K+1),...
%                    ZEV(2*K+1)+0.084*fugou,NONC,fid99,fid100); % OK                         
%                    [GW_3,GR_3,GZ_3] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.126*fugou,REV(2*K),ZEV(2*K)+0.126*fugou,REV(2*K+1),...
%                    ZEV(2*K+1)+0.126*fugou,NONC,fid99,fid100); % OK                         
%                    [GW_4,GR_4,GZ_4] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1)+0.168*fugou,REV(2*K),ZEV(2*K)+0.168*fugou,REV(2*K+1),...
%                    ZEV(2*K+1)+0.168*fugou,NONC,fid99,fid100); % OK
%                    GW = GW + GW_1 + GW_2 + GW_3 + GW_4;
%                    GR = GR + GR_1 + GR_2 + GR_3 + GR_4;
%                    GZ = GZ + GZ_1 + GZ_2 + GZ_3 + GZ_4;
%                end
%                if (abs(ZEV(2*K)) < 0.001 )
%                    GW = GW*0;
%                    GR = GR*0;
%                    GZ = GZ*0;
%                end
               for JJ = 1:3
                   KJ2 = 3*(K-1)+JJ;
                   % �Q�d����^�̉��Z
                   PSIS  = PSIS  + GW(JJ)*FFOUT(sum(NCCN)*2+KJ2)*AMYU0;
	               PSISA = PSISA + GR(JJ)*FFOUT(sum(NCCN)*2+KJ2)*AMYU0;
	               PSISB = PSISB + GZ(JJ)*FFOUT(sum(NCCN)*2+KJ2)*AMYU0;
               end
           end      


%            K = 1:KNE;
%            [GW(1:3,K),GR(1:3,K),GZ(1:3,K)] =  EXTINDC(AS,BS,REV(2*K-1),ZEV(2*K-1),REV(2*K),ZEV(2*K),...
%            REV(2*K+1),ZEV(2*K+1),NONC,fid99,fid100); % OK
%            JJ = 1:3;
%            KJ2(K,1) = 3*(K-1)+1;
%            KJ2(K,2) = 3*(K-1)+2;
%            KJ2(K,3) = 3*(K-1)+3;
%            PSIS  = PSIS + sum(sum(GW(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
% 	       PSISA = PSISA+ sum(sum(GR(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
% 	       PSISB = PSISB+ sum(sum(GZ(JJ,K).*FFOUT(NCCN*2+KJ2(K,JJ))*AMYU0));
%               
                

       end
   else
   end
  %  
  %    ���艻��̉Q�d�������
   if (KSE > 0)  
       for K = 1:KSE
           [GW,GR,GZ] =  EXTINDC(AS,BS,RES(2*K-1),ZES(2*K-1),RES(2*K),ZES(2*K),...
           RES(2*K+1),ZES(2*K+1),NONC,fid99,fid100); % OK
           for JJ = 1:3
	           KJ2 = 2*K-2+JJ;
	           PSIS  = PSIS  + GW(JJ)*FFOUT(sum(NCCN)*2+KNN+KJ2)*AMYU0;
         	   PSISA = PSISA + GR(JJ)*FFOUT(sum(NCCN)*2+KNN+KJ2)*AMYU0;
	           PSISB = PSISB + GZ(JJ)*FFOUT(sum(NCCN)*2+KNN+KJ2)*AMYU0;
           end
       end    
   else
   end
%%%c        write(IPR,*) PSIS,PSISA,PSISB
% ??????????????????????????????????????????????????????????????????
% PSIS�F�C�ӓ_�ł̃��i�O���R�C����^�{�Q�d����^�{CCS��^�j
    for III = 1:CCS
       for K = 1:NE(III)
           [HW,GW,GR,GZ,HR,HZ] = INTEGS(AS,BS,RCCS(III,2*K-1),ZCCS(III,2*K-1),RCCS(III,2*K),...
           ZCCS(III,2*K),RCCS(III,2*K+1),ZCCS(III,2*K+1)); % OK
           for JJ = 1:3
               KJ2 = 3*(K-1)+JJ + 3*sum(NE(1:III-1));
               %�����
               PSIS = PSIS + DFI(KJ2)*GW(JJ)-FI(KJ2)*HW(JJ);
               PSISA= PSISA+ DFI(KJ2)*GR(JJ)-FI(KJ2)*HR(JJ);
               PSISB= PSISB+ DFI(KJ2)*GZ(JJ)-FI(KJ2)*HZ(JJ);
           end
       end
   end
%     K = 1:NE
%     [HW,GW,GR,GZ,HR,HZ] = INTEGS(AS,BS,RCCS(2*K-1),ZCCS(2*K-1),RCCS(2*K),...
%     ZCCS(2*K),RCCS(2*K+1),ZCCS(2*K+1)); % OK
%         for JJ = 1:3
%             KJ2 = 3*(K-1)+JJ;
%             PSIS = PSIS + DFI(KJ2)*GW(JJ)-FI(KJ2)*HW(JJ);
%             PSISA= PSISA+ DFI(KJ2)*GR(JJ)-FI(KJ2)*HR(JJ);
%             PSISB= PSISB+ DFI(KJ2)*GZ(JJ)-FI(KJ2)*HZ(JJ);
%         end
%     end
end
%
%% EDDYP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [n50] = EDDYP(FFOUT,n50,KNE,KNN,NCCN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,REVN,ZEVN)
DISI = zeros(1,594);
EDI = zeros(1,594);
%
% ###############################################################
%  �^��e��ɉ������Q�d�����z�̏o��   (*** #4 ***)
% ###############################################################
      fid40 = fopen('EddyCurrentProfile.txt','w');
      if (NONC <= 0)%   GOTO 507
          DIS = 0.0;
          fprintf(fid40,'%d %d\r\n',DIS,FFOUT(sum(NCCN)*2+1));
          for I = 2:KNN
              SEG = (REV(I)-REV(I-1))^2+(ZEV(I)-ZEV(I-1))^2;
              SEG = sqrt(SEG);
              DIS = DIS+SEG;
              fprintf(fid40,'%d %d\r\n', DIS,FFOUT(sum(NCCN)*2+I));
          end
          SEG = (REV(1)-REV(KNN))^2+(ZEV(1)-ZEV(KNN))^2;
          SEG = sqrt(SEG);
          DIS = DIS+SEG;
          fprintf(fid40,'%d %d\r\n',DIS,FFOUT(sum(NCCN)*2+1));
      else
      % Not yet modified
          II = 0;
          DIS = 0.0D0;
          TOTR = .0D0;
          for K = 1:KNE
              II = II+1;
              SEG1 = (REV(2*K-1)-REVN(3*K-2))^2+(ZEV(2*K-1)-ZEVN(3*K-2))^2;
              DIS = DIS+sqrt(SEG1);
              TOTR = TOTR+FFOUT(sum(NCCN)*2+II)*sqrt(SEG1);
              fprintf(fid40,'%d %d\r\n', DIS,FFOUT(sum(NCCN)*2+II));
              II = II+1;
              SEG2 = (REVN(3*K-2)-REVN(3*K-1))^2+(ZEVN(3*K-2)-ZEVN(3*K-1))^2;
              DIS = DIS+sqrt(SEG2);
              TOTR = TOTR+FFOUT(sum(NCCN)*2+II)*sqrt(SEG2);
              fprintf(fid40,'%d %d\r\n',DIS,FFOUT(sum(NCCN)*2+II));
              II = II+1;
              SEG3 = (REVN(3*K-1)-REVN(3*K))^2+(ZEVN(3*K-1)-ZEVN(3*K))^2;
              DIS = DIS+sqrt(SEG3);
              TOTR = TOTR+FFOUT(sum(NCCN)*2+II)*sqrt(SEG3);
              fprintf(fid40,'%d %d\r\n',DIS,FFOUT(sum(NCCN)*2+II));
              SEG4 = (REVN(3*K)-REV(2*K+1))^2+(ZEVN(3*K)-ZEV(2*K+1))^2;
              DIS = DIS+sqrt(SEG4);
          end
              fprintf('%s %d\r\n','Total reconstructed eddy current =',TOTR);
%
              fid42 = fopen('jeddy0803.txt','r');
              fid43 = fopen('Reference_EddyCurrentProfile.txt','w');
              JMAX = 594;
              TOTI = 0.0;
              DISMAE = 0.0D0;
%              temp_m = zeros(JMAX);
              temp = textscan(fid42,'%s','Delimiter','\n');
              for J = 1:JMAX
                 temp_m = strsplit(temp{1}{J});
                 DISI(J) = str2double(temp_m(1));
                 EDI (J) = str2double(temp_m(2));
                 TOTI=TOTI+(DISI(J)-DISMAE)*EDI(J);
                 DISMAE=DISI(J);
              end
              fprintf('%s %d\r\n','Total refrerence eddy current =',TOTI);
              FACT = TOTR/TOTI;
              J = 1:JMAX;
              EDI(J) = EDI(J)*FACT;
              fprintf(fid43,'%d %d\r\n',horzcat(DISI(J)',EDI(J)'));
%
              fid41 = fopen('EddyCurrent_DiscontinuousProfile.txt','w');
              NSEG = 100;
              ASEG = NSEG;
              DEL = 2.0D0/ASEG;
              DIS = 0.0D0;
              RMAE = REV(1);
              ZMAE = ZEV(1);
              for K = 1:KNE
                  for I = 1:NSEG+1
                      GII = -1.0D0+DEL*(I-1);
                      %-----------------------------------------------------------------------
                      %C ���}�֐��Ă̌v�Z
                      % ZETA1:��1=(3/4)��((3/2)��-1) 
                      % ZETA2:��2=(1-(3/2)��)*(1+(3/2)��)
                      % ZETA3:��3=(3/4)��((3/2)��+1) 
                      %
	                  ZETA1 = 0.75D0*GII*(1.5D0*GII-1.0D0); 
	                  ZETA2 = (1.0D0-1.5D0*GII)*(1.0D0+1.5D0*GII);
	                  ZETA3 = 0.75D0*GII*(1.5D0*GII+1.0D0); 
	                  R = REVN(3*K-2)*ZETA1+REVN(3*K-1)*ZETA2+REVN(3*K)*ZETA3;
	                  Z = ZEVN(3*K-2)*ZETA1+ZEVN(3*K-1)*ZETA2+ZEVN(3*K)*ZETA3;
	                  EDDY = FFOUT(sum(NCCN)*2+3*K-2)*ZETA1+FFOUT(sum(NCCN)*2+3*K-1)*ZETA2+FFOUT(sum(NCCN)*2+3*K)*ZETA3;
	                  DIS = DIS+sqrt((R-RMAE)^2+(Z-ZMAE)^2); 
                      fprintf(fid41,'%d %d\r\n',DIS,EDDY);
                	  RMAE = R;
	                  ZMAE = Z;
                  end
                  fprintf(fid41,'%s\r\n', '==');
              end
      end
end
%

%% SADDLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [statR,endR,AAA1,BBB1,PPPSI1,LOCAT,CCSNEWTON] = SADDLE(CR,CZ,PSI,Z0,RR,...
         FFOUT,n100,n50,RCCS,ZCCS,RC,ZC,ECI,KCMX,RCCN,ZCCN,NCCN,...
         KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,GETA,WAHAHA,MXCCS,...
         RMYU0,NE,NCCS,fid99,fid100,CCS,AUTO,AUTOINPUT)
NCX = 200;
NCY = 200;
CXK = zeros(NCX,NCY);
CYK = zeros(NCX,NCY);
SOLK = zeros(NCX,NCY);
FI = zeros(1,MXCCS);
DFI = zeros(1,MXCCS);%  ! BOUNDARY CONDITION FI:��, DFI:d��/dn
ECIGRP = 0; % ushiki
AAA1 = 0;% ushiki
BBB1 = 0;% ushiki
statR = 0;% ushiki
endR =  0;% ushiki
CCSNEWTON =1;
%
      fprintf('%s\r\n','##');
      fprintf('%s\r\n','## �ŊO�k���C�ʂ܂��̓v���Y�}�d�����E�̒T�� ##');
      fprintf('%s\r\n','�����֐��̓��������z������B');
      if (AUTO == 0)
          prompt = '���C�ʏ�̓_���w��? / �w�_(Bp=0)���T�[�`? =(1/0)\n';
          LOCAT = input(prompt);
      else
          LOCAT = AUTOINPUT(37);
      end
      if (LOCAT > 0)
          if (CCS > 1) 
              prompt = '���Ԗڂ�CCS���S����NEWTON�@���s��? \n';
              CCSNEWTON = input(prompt);
          end
          if (AUTO == 0)
              fprintf('���C�ʏ�̓_�̍��W(R,Z)�́H\n');
              prompt = 'R = ';
              AAA1 = input(prompt);
              prompt = 'Z = ';
              BBB1 = input(prompt);
%               prompt = '�ŊO�k���C��(LCMS)? / �v���Y�}�d�����E�H =(0/1)\n';
%           LCBP = input(prompt);
          else
              AAA1 = AUTOINPUT(38+ (CCSNEWTON-1)*2);
              BBB1 = AUTOINPUT(39+ (CCSNEWTON-1)*2);              
          end
      else
          if (AUTO == 0)
              prompt = 'X�_�T�[�`�̍����ʒu(r[m])\n';
	          statR = input(prompt);
              prompt = 'X�_�T�[�`�̉E���ʒu(r[m])\n';
              endR = input(prompt);
          else
	          statR = AUTOINPUT(42);
              endR  = AUTOINPUT(43);              
          end
          fprintf(WAHAHA,'%s\r\n','X�_�T�[�`');
          %!******** SERCH NULL POINT ***************************************
          IR = 1; % 1000
          IZ = 0;
	      CRMAE = CR(1);
          for I = 1:NINT
	          DCR = abs(CR(I)-CRMAE);
%              if (DCR < 1.0D-10) %GOTO 49
%              else
	              IR = IR+1.*(DCR >= 1.0D-10);
	              IZ = IZ.*(DCR < 1.0D-10);
%              end
              IZ = IZ+1;
	          CXK(IR,IZ) = CR(I);
	          CYK(IR,IZ) = CZ(I);
	          SOLK(IR,IZ) = PSI(I);
	          CRMAE = CR(I);
          end
          IRMX = IR; 
    	  IZMX = IZ;
      end
      RCCS(NCCS+1) = RCCS(1); % 1999
      ZCCS(NCCS+1) = ZCCS(1);
%%%%%%C	DO I=1,NCCS
      DFI(1:NCCN) = FFOUT(1:NCCN);
%%%%%%%C      DO I=1,NCCS
%%%%%%%%C        FI(I)=FFOUT(I+NCCS)
      FI(1:NCCN) = FFOUT(1 + NCCN:NCCN + NCCN);
%%
      if (LOCAT > 0)
      else
          BELOW = Z0-RR;
	      BPMIN = 1.0D20;
          for IR=1:IRMX
              for IZ = 1:IZMX
                  A1 = CXK(IR,IZ);
                  B1 = CYK(IR,IZ);
                  if(CYK(IR,IZ) > BELOW)
                      continue
                  end
                  if or(CXK(IR,IZ) < statR, CXK(IR,IZ) > endR)
                      continue
                  end
%
% *******************************************************************
                  [PSI0,PSIA,PSIB] = QINTER(A1,B1,GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
                  n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,...
                  REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
% ********************************************************************CC 
                  BP2 = (PSIA/A1)^2+(PSIB/A1)^2;
                  EEE = PSIA/abs(PSIA)+PSIB/abs(PSIB);
                  EEE = abs(EEE);
%?????????????????????????????????????????      
                  %cc	IADD=1 
	              IADD = 0; 
%?????????????????????????????????????????      
                  if (IADD > 0)
                      BP2 = BP2+EEE;
                  end
%+++
                  if(BP2 > BPMIN)
                      continue
                  end
	              BPMIN = BP2;
	              IRRR = IR;
	              IZZZ = IZ;
              end%100
          end%110
          fprintf(WAHAHA,'%d X-point=%d Bp2=%d Flux=%d\r\n', CXK(IRRR,IZZZ),...
          CYK(IRRR,IZZZ),BPMIN,SOLK(IRRR,IZZZ));
%
% insert��
%
          RMIN = CXK(IRRR-1,IZZZ);
	      RMAX = CXK(IRRR+1,IZZZ);
	      ZMIN = CYK(IRRR,IZZZ-1);
	      ZMAX = CYK(IRRR,IZZZ+1);
%         220
          ISECT = 10;
	      CSECT = ISECT;
	      DELR = (RMAX-RMIN)/CSECT;
	      DELZ = (ZMAX-ZMIN)/CSECT;
	      A1 = RMIN;
          for IR = 1:ISECT-1
	          A1 = A1+DELR;
	          B1 = ZMIN;
              for IZ = 1:ISECT-1
	              B1 = B1+DELZ;
%
% *******************************************************************
                  [PSI0,PSIA,PSIB] = QINTER(A1,B1,GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
                  n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,...
                  REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
%********************************************************************CC 
                  BP2 = (PSIA/A1)^2+(PSIB/A1)^2;
                  EEE = PSIA/abs(PSIA)+PSIB/abs(PSIB);
	              EEE = abs(EEE);
                  if (IADD > 0)
                      BP2 = BP2+EEE;
                  end
                  if (BP2 > BPMIN)
                      continue
                  end
	              BPMIN = BP2;
	              AAA1 = A1;
	              BBB1 = B1;
               
              end
% insert��
          
          end
      end
% *******************************************************************
      [PSI1,PSIA,PSIB] = QINTER(AAA1,BBB1,GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
      n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,...
      RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
% ********************************************************************CC 
%
      if (LOCAT <= 0)        
          fprintf(WAHAHA,'%d %d %d %d\r\n',AAA1,BBB1,BPMIN,PSI1);
          RMIN = AAA1-DELR;
          RMAX = AAA1+DELR;
	      ZMIN = BBB1-DELZ;
	      ZMAX = BBB1+DELZ;
          while ( DELR*DELZ > 1.0E-10) % goto 220
              ISECT = 10;
	          CSECT = ISECT;
	          DELR = (RMAX-RMIN)/CSECT;
	          DELZ = (ZMAX-ZMIN)/CSECT;
	          A1 = RMIN;
              for IR = 1:ISECT-1
	              A1 = A1+DELR;
	              B1 = ZMIN;
                  for IZ = 1:ISECT-1
	                  B1 = B1+DELZ;
%
% *******************************************************************
                      [PSI0,PSIA,PSIB] = QINTER(A1,B1,GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
                      n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,...
                      REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
%********************************************************************CC 
                      BP2 = (PSIA/A1)^2+(PSIB/A1)^2;
                      EEE = PSIA/abs(PSIA)+PSIB/abs(PSIB);
	                  EEE = abs(EEE);
                      BP2 = BP2 + EEE.*(IADD > 0);
                      if (BP2 > BPMIN)
                          continue
                      end
	                  BPMIN = BP2;
	                  AAA1 = A1;
	                  BBB1 = B1;
                  end
              end
      
% *******************************************************************
              [PSI1,PSIA,PSIB] = QINTER(AAA1,BBB1,GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
              n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,...
              RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
% ********************************************************************CC 
%  
              fprintf(WAHAHA,'%d %d %d %d\r\n',AAA1,BBB1,BPMIN,PSI1);
              RMIN = AAA1-DELR;
              RMAX = AAA1+DELR;
	          ZMIN = BBB1-DELZ;
	          ZMAX = BBB1+DELZ;
          end
      end
      fprintf('NULL POINT %d %d  PSI = %d\r\n',AAA1,BBB1,PSI1);
      PPPSI1 = PSI1;
end
%
%% NEWTON!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function[TOTMOD,TRUTOT,EPSAX,DEFR,DEFZ,C0R,C0Z,I1I0,ELNG] = NEWTON(IUTST,ICONT,RADCCS,RCCS,ZCCS,...
FFOUT,n100,n50,RR0,ZZ0,RC,ZC,ECIGRP,ECI,KCMX,A1,B1,PSI1,TOTMOD,I1I0,RCCN,ZCCN,NCCN,...
SOU,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,GETA,LOCAT,NCCS,RMYU0,fid99,fid100,fid50,fid51,WAHAHA,MXCCS,CCS,NE,CCSNEWTON,AUTO,AUTOINPUT)
FI  = zeros(1,MXCCS);
DFI = zeros(1,MXCCS); % ! BOUNDARY CONDITION FI:��, DFI:d��/dn
SP = zeros(1,1300);
ABCR = zeros(1,1300);
ABCZ = zeros(1,1300);%,
DEFR = zeros(1,1300);
DEFZ = zeros(1,1300);
DNPSI = zeros(1,1300);
%
% NE = NCCS/2; ushiki!
%%
for III = 1:CCS
    RCCS(III,NCCS+1)=RCCS(III,1);
    ZCCS(III,NCCS+1)=ZCCS(III,1);
end
%RCCS(NCCS+1) = RCCS(1);
%ZCCS(NCCS+1) = ZCCS(1);
DFI(1:sum(NCCN)) = FFOUT(1:sum(NCCN));
FI(1:sum(NCCN)) = FFOUT(1 + sum(NCCN):sum(NCCN) + sum(NCCN));
%DFI(1:NCCN) = FFOUT(1:NCCN);
%FI(1:NCCN) = FFOUT(1+NCCN:NCCN+NCCN);
%!*********************************************** 
%!**********  NEWTON'S METHOD  ****************** OK
%!***********************************************
fid18 =	fopen('BP_VECTOR.TXT','w');
fid29 = fopen('DISTANCE.TXT','w');
fid30 = fopen('INVERSE_input.TXT','w');
	frewind(fid18);%
	frewind(fid29);
	frewind(fid30);%
%
    fprintf(WAHAHA, '%s\r\n','Newton�@�ɂ��ŊO�k���C��');
%cocococococo
    NB = 300;%  !�ŊO�k���C�ʂ̓_�̐�  
%cocococococo
    RNB = NB;%
    AS = A1;
    BS = B1;
    PSIS = PSI1;
    fprintf(WAHAHA,'%d %d %d\r\n', AS,BS,PSIS);
% ********************************************************************
    [PSIS,PSISA,PSISB] = QINTER(AS,BS,GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,...
    n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,...
    Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
% ********************************************************************
%
% CCS����̏ꍇ�ɂ�CCS���S���v�Z��������悤�ɉ��ǂ���K�v����
if (CCS > 1)
%     prompt = '���Ԗڂ�CCS���S����NEWTON�@���s��? \n';
%     CCSNEWTON = input(prompt);
    RR0 = RR0(CCSNEWTON);
    ZZ0 = ZZ0(CCSNEWTON);
else
    CCSNEWTON = 1;
end
if and(LOCAT > 0, CCSNEWTON == 2)
	frewind(fid51);%!  (51,FILE='BOUNDARY_P.TXT')
else
	frewind(fid50);%!  (50,FILE='BOUNDARY.TXT')
end
    R0 = sqrt((AS-RR0)*(AS-RR0)+(BS-ZZ0)*(BS-ZZ0));   
    FAI0 = -pi+asin((ZZ0-BS)/R0); %      !+ �� - ?
%
    AR = RR0+R0*cos(FAI0);
    BZ = ZZ0+R0*sin(FAI0);
	fprintf(WAHAHA,'AR=%d RZ=%d TET=%d\r\n', AR,BZ,FAI0);
%%
    DFAI = 2.0D0*pi/RNB;
    NCOUNT = 1;
    fprintf(WAHAHA, '#  R0 = %d',R0);
    fprintf(WAHAHA,'TET = %d Rad = %d R/Z = %d %d PSI = %d\r\n',FAI0,R0,NCOUNT,AS,BS,PSIS);
    IWR = 1;
	RWR(IWR) = AS;
	ZWR(IWR) = BS;
	DNPSI(IWR) = sqrt(PSISA*PSISA+PSISB*PSISB);
    fprintf(WAHAHA, '%d %d %d %d\r\n', IWR,PSISA,PSISB,DNPSI(IWR));
%     if and(LOCAT > 0, LCBP > 0)
%         fprintf(fid51,'%d %d\r\n', AS,BS);%  ! ##2
%     else
%         fprintf(fid50,'%d %d\r\n', AS,BS);%  ! ##2
%     end
    FAI1 = FAI0+DFAI;
 
    
    %
    [IWR,NCOUNT,R3HIT,R3OLD,RWR,ZWR,R0,AMIN,AMAX,DNPSI,R1,FAI1,F1,F2] = SERCH_LCMS...
    (GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN...
    ,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,R0,RR0,...
    FAI1,IWR,NCOUNT,ZZ0,fid99,fid100,RMYU0,NE,RADCCS,PSIS,WAHAHA,LOCAT,fid50,fid51,...
    DFAI,DNPSI,RWR,ZWR,0,0,CCS,CCSNEWTON);
%
    while (NCOUNT ~= NB) 
 
        [IWR,NCOUNT,R3HIT,R3OLD,RWR,ZWR,R0,AMIN,AMAX,DNPSI,R1,FAI1,F1,F2] = SERCH_LCMS(GETA,...
        RCCS,ZCCS,FFOUT,FI,DFI,n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,...
        KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,R0,RR0,FAI1,IWR,NCOUNT,ZZ0,...
        fid99,fid100,RMYU0,NE,RADCCS,PSIS,WAHAHA,LOCAT,fid50,fid51,DFAI,...
        DNPSI,RWR,ZWR,F1,F2,CCS,CCSNEWTON);
%    
    end
    fprintf(WAHAHA, '%d %d %d %d %d %d\r\n',FAI1,R1,NCOUNT,AS,BS,PSIS);
%+++
    IWRMX = IWR;
    IWR = IWR+1;
    RWR(IWR) = AS;
	ZWR(IWR) = BS;
	DNPSI(IWR) = sqrt(PSISA*PSISA+PSISB*PSISB);
    fprintf(WAHAHA,'%d %d %d %d\r\n', IWR,PSISA,PSISB,DNPSI(IWR));
%+++  
    if and(LOCAT > 0, CCSNEWTON == 2)
        fprintf(fid51,'%d %d\r\n', RWR(2),ZWR(2)); % !  ##4
    else
        fprintf(fid50,'%d %d\r\n', RWR(2),ZWR(2)); % !  ##4
    end
    fprintf('%s\r\n','OK')
    RWR(IWR) = [];
	ZWR(IWR) = [];
	DNPSI(IWR) = [];
    RWR(1) = [];
	ZWR(1) = [];
	DNPSI(1) = [];
    IWRMX = IWRMX -2;
    IWR = IWR -2;
%
    DNP = zeros(1,IWRMX);
    IWR = 1:IWRMX;
    DNP(IWR) = (DNPSI(IWR)+DNPSI(IWR+1))./2.0D0;
    fprintf(WAHAHA, '%d %d %d DNPSI = %d DNP = %d', horzcat(IWR',RWR(IWR)',...
                                         ZWR(IWR)',DNPSI(IWR)',DNP(IWR)'));
%
    IWR = 1:IWRMX;
    fprintf(WAHAHA,'%d %d %d\r\n',horzcat(IWR',RWR(IWR)',ZWR(IWR)'));
    fprintf(fid30,'%d %d %d\r\n',horzcat(IWR',RWR(IWR)',ZWR(IWR)'));
%
    II3 = 3;
    PP0 = PSIS;
      for IWR = 1:IWRMX
	      DNP1 = -(DNPSI(IWR)+DNPSI(IWR+1))/2.0D0;
          fprintf(WAHAHA,'%d %d %d %d %d\r\n', IWR,IWR,II3,PP0,DNP1);
      end
      RWR(IWRMX+1) = RWR(1);
      ZWR(IWRMX+1) = ZWR(1);
      R0 = RWR(1);
	  Z0 = ZWR(1);
	  DIST = 0.0D0;
      for  IWR = 1:IWRMX+1
	      DELR = RWR(IWR)-R0;
	      DELZ = ZWR(IWR)-Z0;
	      DIST = DIST+sqrt(DELR*DELR+DELZ*DELZ);
	      BP = DNPSI(IWR)/RWR(IWR);
	      fprintf(fid29,'%d %d\r\n', DIST,BP);
          R0 = RWR(IWR);
	      Z0 = ZWR(IWR);
      end
%
      RWR(IWRMX+1) = RWR(1);
      ZWR(IWRMX+1) = ZWR(1);
	  TOT = 0.0D0;
%
      fprintf('%s %d\r\n','TOTMOD =',TOTMOD);
%
      for IWR = 1:IWRMX
          A3 = (RWR(IWR)+RWR(IWR+1))/2.0d0;
          B3 = (ZWR(IWR)+ZWR(IWR+1))/2.0d0;
          DELR = RWR(IWR+1)-RWR(IWR);
          DELZ = ZWR(IWR+1)-ZWR(IWR);
	      DIST = sqrt(DELR*DELR+DELZ*DELZ);
% *******************************************************************
          [PSI3,PSI3A,PSI3B] = QINTER(A3,B3,GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100...
    ,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,...
    ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
% ********************************************************************CC
          II3 = 3;
	      PP0 = PSIS;
	      DNP1 = -sqrt(PSI3A*PSI3A+PSI3B*PSI3B);
%
          AR = RWR(IWR+1)-RWR(IWR);
          AZ = ZWR(IWR+1)-ZWR(IWR);
	      AL = sqrt(AR*AR+AZ*AZ);
	      DNR = AZ/AL;
	      DNZ = -AR/AL;
	      DNP2 = -(DNR*PSI3A+DNZ*PSI3B);
%
          DNP2 = DNP2*TOTMOD;
%
          fprintf(WAHAHA,'%d %d %d %d %d\r\n', IWR,IWR,II3,PP0,DNP2);
          fprintf(fid30, '%d %d %d %d %d\r\n',IWR,IWR,II3,PP0,DNP2);
%
    	  EXPND=3.0D0;
	      PSI3AR=-PSI3B/A3;
	      PSI3BZ=PSI3A/A3;
	      A33=A3+EXPND*PSI3AR;
	      B33=B3+EXPND*PSI3BZ;
%	
%******
%
          TETP=atan2(PSI3BZ,PSI3AR);
          PSN0=PSI3AR*cos(TETP)+PSI3BZ*sin(TETP);
          PSN=abs(PSN0);%   ! ����
%     ���ǁAPSN=DSQRT(PSI3AR*PSI3AR+PSI3BZ*PSI3BZ)�ƌ��ʂ͕ς��Ȃ��B      
%******
%
	      TOT=TOT+PSN*DIST;
          fprintf(fid18,'%d %d %d %d\r\n', A3,B3,A33,B33);

      end
      AMYU0 = 4.0D0*pi*1.0D-1;%! NAMUAMUdabutsu  #3
	  TOT = TOT/AMYU0;
      fprintf('%s   %d       %s\r\n','Total current =',TOT, '(MA)');
      if (ICONT > 1)
          TOTMOD = 1.0D0.*(I1I0 <= 0) + TOTMODif.*(I1I0 > 0);
      else
          if (AUTO == 0)
    	      prompt =  'Modify the normal derivative output? (Yes/No)=(1/0)\n';
              I1I0 = input(prompt);
          else
              I1I0 = AUTOINPUT(44);
          end
          if (I1I0 <= 0)
              TOTMOD = 1.0D0;
          elseif (ICONT > 1)
              TOTMOD = TRUTOT/TOT;
              fprintf('%s   %d       %s\r\n','Total current now modified as TRUTOT=',TRUTOT,'(MA)');
          else
              if (AUTO == 0)
	              prompt = 'Input the true value of total current.\n';
	              TRUTOT = input(prompt);
              else
                  TRUTOT = AUTOINPUT(45 + (CCSNEWTON-1)*1);
              end
              TOTMOD = TRUTOT/TOT;
              fprintf('%s   %d       %s\r\n','Total current now modified as TRUTOT=',TRUTOT,'(MA)');
% *******************************************************************
%   CCS���S�̏C��
          end
      end
      CTOP = -1.0D10;
	  CBTM =  1.0D10;
	  CLFT =  1.0D10;
	  CRIT = -1.0D10;
      for  IWR = 1:IWRMX
          CTOP = ZWR(IWR).*(ZWR(IWR) > CTOP) + CTOP.*(ZWR(IWR) <= CTOP);
          CBTM=ZWR(IWR).*(ZWR(IWR) < CBTM) + CBTM.*(ZWR(IWR) >= CBTM);
          CLFT=RWR(IWR).*(RWR(IWR) < CLFT) + CLFT.*(RWR(IWR) >= CLFT);
          CRIT=RWR(IWR).*(RWR(IWR) > CRIT) + CRIT.*(RWR(IWR) <= CRIT);
      end
%  ���E�ϕ��ɂ��d�S�̌���
%  
    RWR(IWRMX+1) = RWR(1);
    ZWR(IWRMX+1) = ZWR(1);
    TOTA = 0.0D0;
    TOTR = 0.0D0;
    TOTZ = 0.0D0;
    WAY = 0.0D0;
    for IWR = 1:IWRMX
        RM = (RWR(IWR+1)+RWR(IWR))*0.5D0;
        ZM = (ZWR(IWR+1)+ZWR(IWR))*0.5D0;
        ER = RWR(IWR+1)-RWR(IWR);
        EZ = ZWR(IWR+1)-ZWR(IWR);
        SL = sqrt(ER*ER+EZ*EZ);
        RNOR = EZ/SL;
        ZNOR = -ER/SL;
        EET = (RNOR*RM+ZNOR*ZM)*0.5D0;
        EER = RNOR;
        EEZ = ZNOR;
        TOTA = TOTA+EET*SL;
        PHI = (RM*RM+ZM*ZM)*0.25D0;
        TOTR = TOTR+(RM*EET-PHI*EER)*SL;
        TOTZ = TOTZ+(ZM*EET-PHI*EEZ)*SL;
        WAY = WAY+SL;
    end
        C0R = TOTR/TOTA; % !�@�d�S��r���W
        C0Z = TOTZ/TOTA ; %!�@�d�S��z���W
        fprintf('%s %d %d\r\n', 'NEW CCS CENTER',C0R,C0Z);
        fprintf('%s %d\r\n', '���̂�', WAY);
%
        WAYX = WAY/NCCS(CCSNEWTON);% ! ���̂�/NCCS
        SP(1:NCCS(CCSNEWTON)) = WAYX*(1:NCCS(CCSNEWTON));
        for I = 1:NCCS(CCSNEWTON)-1
            WAY=0.0D0;
            for IWR = 1:IWRMX
                ER = RWR(IWR+1)-RWR(IWR);
                EZ = ZWR(IWR+1)-ZWR(IWR);
                SL = sqrt(ER*ER+EZ*EZ);
                WAY = WAY+SL;
                if(WAY >= SP(I))
                    ABCR(I+1) = RWR(IWR);% ! �k���O��LCMS��r���W�iCCS���b�V���_�����ӏ��j
                    ABCZ(I+1) = ZWR(IWR); %! �k���O��LCMS��z���W�iCCS���b�V���_�����ӏ��j
                    break
                end
            end
        end
        ABCR(1) = RWR(1); %! �k���O��LCMS��r���W�iX�_�ʒu�j
        ABCZ(1) = ZWR(1) ;%! �k���O��LCMS��z���W�iX�_�ʒu�j
%
%   LCMS�ɑ�����CCS�`��̌���i�k���� SOU�j
%   DEFR(I),DEFZ(I) = �VCCS�̃��b�V���_���W (I=1,NCCS)
%        
	ELNG = (CTOP-CBTM)/(CRIT-CLFT);
    fprintf('%s %d %s %d\r\n', 'CTOP=',CTOP,'  CBTM=',CBTM);
	fprintf('%s %d %s %d\r\n','CLFT=',CLFT,'  CRIT=',CRIT);
	fprintf('%s %d %d\r\n','Co-Center=',C0R,C0Z);
	fprintf('%s %d\r\n','TateYoko-Ratio=',ELNG);
	EPSAX = sqrt((C0R-RR0)^2+(C0Z-ZZ0)^2);
	fprintf('%s %d\r\n','Deviation OF CCS-Center',EPSAX);
	DEFR(1:NCCS(CCSNEWTON)) = C0R+(ABCR(1:NCCS(CCSNEWTON))-C0R)/SOU(CCSNEWTON);
	DEFZ(1:NCCS(CCSNEWTON)) = C0Z+(ABCZ(1:NCCS(CCSNEWTON))-C0Z)/SOU(CCSNEWTON);
	ABCR(NCCS(CCSNEWTON)+1) = ABCR(1);
	ABCZ(NCCS(CCSNEWTON)+1) = ABCZ(1);
	DEFR(NCCS(CCSNEWTON)+1) = DEFR(1);
	DEFZ(NCCS(CCSNEWTON)+1) = DEFZ(1);
end
%
%%   SERCH_LCMS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
function [IWR,NCOUNT,R3HIT,R3OLD,RWR,ZWR,R0,AMIN,AMAX,DNPSI,R1,FAI1,F1,F2] = SERCH_LCMS...
         (GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,...
         REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,R0,RR0,FAI1,IWR,NCOUNT,ZZ0,fid99,fid100,RMYU0,...
         NE,RADCCS,PSIS,WAHAHA,LOCAT,fid50,fid51,DFAI,DNPSI,RWR,ZWR,F1,F2,CCS,CCSNEWTON)
R3X = zeros(1,21);
FX3X = zeros(1,21);
FX1 = 0;% ??????????????? ushiki

   %C
   %CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   %C    Newton�@�ɂ��LCMS��̓_���T�[�`  (�� START)
   %CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   %C   
   %C�����ł��B���ڂ��Ă��������B�i���j         
    RLIMT = 0.10815D0;
   %C�����ł��B���ڂ��Ă��������B�i���j  
    RMIN = R0*0.98;
    RMAX = R0*1.02;
    AMIN = RR0+RMIN*cos(FAI1);
    AMAX = RR0+RMAX*cos(FAI1);
    LXEED = 0;
    LXEED = 1.*(AMAX < RLIMT) + LXEED.*(AMAX >= RLIMT);
    LXEED = 2.*(AMIN < RLIMT) + LXEED.*(AMIN >= RLIMT);
    RMAX = (RLIMT-RR0)/cos(FAI1).*(LXEED > 0) + RMAX.*(LXEED >= 0);
    RMIN = (RLIMT-0.001D0-RR0)/cos(FAI1).*(LXEED > 0) + RMIN.*(LXEED >= 0);
    AMIN = RR0+RMIN*cos(FAI1).*(LXEED > 0) + AMIN.*(LXEED >= 0);
    AMAX = RR0+RMAX*cos(FAI1).*(LXEED > 0) + AMAX.*(LXEED >= 0);
    NN = 10;
    CN = NN;
    RDEL = (RMAX-RMIN)/CN;
    FXMIN = 1.0D10;
    for JJ = 1:NN+1
        R3X(JJ) = RMIN+RDEL*(JJ-1);
        A3 = RR0+R3X(JJ)*cos(FAI1);
        B3 = ZZ0+R3X(JJ)*sin(FAI1);
        % *******************************************************************
        [PSI3,PSI3A,PSI3B] = QINTER(A3,B3,GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,...
        n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,...
        Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
        % ********************************************************************C
        FX3X(JJ) = PSI3-PSIS;
        if(abs(FX3X(JJ)) < FXMIN)
            FXMIN = abs(FX3X(JJ));
            R3HIT = R3X(JJ);
            JJHIT = JJ;
        end
    end
    if (JJHIT == NN+1)
        R1 = R3X(NN);
        R2 = R3X(NN+1);
        F1 = FX3X(NN);
        F2 = FX3X(NN+1);
        SOLF = F2;
    elseif(JJHIT == 1)
        R1 = R3X(1);
        R2 = R3X(2);
        F1 = FX3X(1);
        F2 = FX3X(2);
        SOLF = F1;
    else
        PM1 =  FX3X(JJHIT-1)*FX3X(JJHIT);
        PM2 = FX3X(JJHIT)*FX3X(JJHIT+1);
        DS1 = abs(FX3X(JJHIT-1))-abs(FX3X(JJHIT));
        DS2 = abs(FX3X(JJHIT))-abs(FX3X(JJHIT+1));
        DS1 = abs(DS1);
        DS2 = abs(DS2);
        if or(PM1 < 0.0D0, DS1 < DS2)
            R1 = R3X(JJHIT-1);
            R2 = R3X(JJHIT);
            F1 = FX3X(JJHIT-1);
            F2 = FX3X(JJHIT);
            SOLF = F2;
        end
        if or(PM2 < 0.0D0, DS2 < DS1)
            R1 = R3X(JJHIT);
            R2 = R3X(JJHIT+1);
            F1 = FX3X(JJHIT);
            F2 = FX3X(JJHIT+1);
            SOLF = F1;
        end
    end
% 333
    LLPRNT = 0;
    R3OLD = 1.0D10;
    for  ICNT = 3:100
        if (ICNT >= 3)
            R3 = R2-F2*(R2-R1)/(F2-F1);
        end
%%%C  Anzen-Souchi
        if or(R3 < RADCCS, R3 > 2.0*R0)
            R3=(R1+R2)/2.0D0;
        end
        A3 = RR0+R3*cos(FAI1);
%%%C�����ł��B���ڂ��Ă��������B�i���j             
        if (A3 < RLIMT)
            LLPRNT = LLPRNT+1;
            RMAX = (RLIMT-RR0)/cos(FAI1); % 19
            RMIN = RMAX-0.05D0;
            NN = 10;
            CN = NN;
            while (1 == 1) %18
                RDIF = abs((RMAX-RMIN)/RMIN);
                if (RDIF < 1.0D-05) 
                    R3 = (RMAX+RMIN)/2.0D0;
                    break
                end
                RDEL = (RMAX-RMIN)/CN;
                FXMIN = 1.0D10;
                for JJ = 1:NN+1
                    R3X(JJ) = RMIN+RDEL*(JJ-1);
                    A3 = RR0+R3X(JJ)*cos(FAI1);
                    B3 = ZZ0+R3X(JJ)*sin(FAI1);
                    % *******************************************************************
                    [PSI3,PSI3A,PSI3B] = QINTER(A3,B3,GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,...
                    n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,...
                    Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
                    % ********************************************************************C
                    FX3X(JJ) = PSI3-PSIS;
                    if (abs(FX3X(JJ)) < FXMIN)
                        FXMIN = abs(FX3X(JJ));
                        R3HIT = R3X(JJ);
                        JJHIT = JJ;
                    end
                end
%  
                if (JJHIT == NN+1)
                    RMIN = R3X(NN);
                    RMAX = R3X(NN+1);
                    F1 = FX3X(NN);
                    F2 = FX3X(NN+1);
                    continue
                end
                if (JJHIT == 1)
                    RMIN = R3X(1);
                    RMAX = R3X(2);
                    F1 = FX3X(1);
                    F2 = FX3X(2);
                    continue
                end
%C      
                PM1 = FX3X(JJHIT-1)*FX3X(JJHIT);
                PM2 = FX3X(JJHIT)*FX3X(JJHIT+1);
                DS1 = abs(FX3X(JJHIT-1))-abs(FX3X(JJHIT));
                DS2 = abs(FX3X(JJHIT))-abs(FX3X(JJHIT+1));
                DS1 = abs(DS1);
                DS2 = abs(DS2);
                if or(PM1 < 0.0D0, DS1 < DS2)
                    RMIN = R3X(JJHIT-1);
                    RMAX = R3X(JJHIT);
                    F1 = FX3X(JJHIT-1);
                    F2 = FX3X(JJHIT);
                    continue
                elseif or(PM2 < 0.0D0, DS2 < DS1)
                    RMIN = R3X(JJHIT);
                    RMAX = R3X(JJHIT+1);
                    F1 = FX3X(JJHIT);
                    F2 = FX3X(JJHIT+1);
                    continue
                else
                    break
                end
            end
            break
        else
            A3 = RR0+R3*cos(FAI1);
            B3 = ZZ0+R3*sin(FAI1);
            %C *******************************************************************
            [PSI3,PSI3A,PSI3B] = QINTER(A3,B3,GETA,RCCS,ZCCS,FFOUT,FI,DFI,n100,...
            n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,...
            Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);
            %C ********************************************************************C
            FX3 = PSI3-PSIS;
            if (LLPRNT == 1)
                fprintf('%s %d %d %d\r\n','ICNT/A3/FX3=',ICNT,A3,FX3);
            end
            if (FX3*FX1<0.0D0)
                R2 = R3;
	            F2 = FX3;
                F1 = F1/2.0D0.*(FX3*SOLF>0.0D0) + F1.*(FX3*SOLF<=0.0D0);
            else
	            R1 = R3;
                F1 = FX3;
%                 if (FX3*SOLF>0.0D0)
% 	                F2 = F2/2.0D0;
%                 end
	            F2 = F2/2.0D0.*(FX3*SOLF>0.0D0) + F2.*(FX3*SOLF<=0.0D0);
            end
	        SOLF = FX3;
            EPSR = 1.0D-8;
            if (abs((R2-R1)/R1) < EPSR) 
                break
            end
            EPSF=1.0D-8;
            if (abs(FX3) < EPSF)
                break
            end
        end
    end
%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
%C    Newton�@�ɂ��LCMS��̓_���T�[�`  (�� END)
%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
%C   
      R1=R3; %21
	  NCOUNT = NCOUNT+1;
	  RRR1 = RR0+R1*cos(FAI1);
	  ZZZ1 = ZZ0+R1*sin(FAI1);
      fprintf(WAHAHA,'%d %d %d %d %d %d\r\n', FAI1,R1,NCOUNT,RRR1,ZZZ1,PSI3);
      IWR = IWR+1;
	  RWR(IWR) = RRR1;
	  ZWR(IWR) = ZZZ1;
	  DNPSI(IWR) = sqrt(PSI3A*PSI3A+PSI3B*PSI3B);
      fprintf(WAHAHA,'%d PSIA = %d PSIB = %d DNPSI = %d\r\n', IWR,PSI3A,PSI3B,DNPSI(IWR));
%+++
      if and(LOCAT > 0, CCSNEWTON == 2)
          fprintf(fid51,'%d %d\r\n', RRR1,ZZZ1);%  ! ##3
      else
          fprintf(fid50,'%d %d\r\n', RRR1,ZZZ1);%  ! ##3
      end
      FAI1=FAI1+DFAI;
      R0=R1;
end
%% MakePlasmaCurrentDensityData!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function MakePlasmaCurrentDensityData(GETA,RCCS,ZCCS,FFOUT,NCCN,NE,CCS,BOUND,MXCCS,IUTST,n100,...
n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,...
RMYU0)
prompt = '��K�����E�v�f��(�ߓ_����1/3)��?�i�����l��40-60�j\n';
NEC = input(prompt);
NMC = 2*NEC;
%%  ���E�ϕ��ɂ��d�S�̌���   
RWR = BOUND(:,1);
ZWR = BOUND(:,2);
IWRMX = numel(BOUND(:,1));
RWR(IWRMX+1) = RWR(1);
ZWR(IWRMX+1) = ZWR(1);
TOTA = 0.0D0;
TOTR = 0.0D0;
TOTZ = 0.0D0;
WAY = 0.0D0;
Mesh = fopen(strcat('UTST_',num2str(IUTST),'/MeshPoint_UTST.txt'),'w');
Node_BC = fopen(strcat('UTST_',num2str(IUTST),'/Node_BoundaryCondition_UTST.txt'),'w');
for IWR = 1:IWRMX
    RM = (RWR(IWR+1)+RWR(IWR))*0.5D0;
    ZM = (ZWR(IWR+1)+ZWR(IWR))*0.5D0;
    ER = RWR(IWR+1)-RWR(IWR);
    EZ = ZWR(IWR+1)-ZWR(IWR);
    SL = sqrt(ER*ER+EZ*EZ);
    RNOR = EZ/SL;
    ZNOR = -ER/SL;
    EET = (RNOR*RM+ZNOR*ZM)*0.5D0;
    EER = RNOR;
    EEZ = ZNOR;
    TOTA = TOTA+EET*SL;
    PHI = (RM*RM+ZM*ZM)*0.25D0;
    TOTR = TOTR+(RM*EET-PHI*EER)*SL;
    TOTZ = TOTZ+(ZM*EET-PHI*EEZ)*SL;
    WAY = WAY+SL;
end
C0R = TOTR/TOTA; % !�@�d�S��r���W
C0Z = TOTZ/TOTA ; %!�@�d�S��z���W
fprintf('%s %d %d\r\n', 'NEW CCS CENTER',C0R,C0Z);
fprintf('%s %d\r\n', '���̂�', WAY);
%
WAYX = WAY/NMC;% ! ���̂�/NCCS
SP(1:NMC) = WAYX*(1:NMC);
XM = zeros(1,NMC+1);
YM = zeros(1,NMC+1);
for I = 1:NMC-1
    WAY=0.0D0;
    for IWR = 1:IWRMX
        ER = RWR(IWR+1)-RWR(IWR);
        EZ = ZWR(IWR+1)-ZWR(IWR);
        SL = sqrt(ER*ER+EZ*EZ);
        WAY = WAY+SL;
        if(WAY >= SP(I))
            XM(I+1) = RWR(IWR);% ! �k���O��LCMS��r���W�iCCS���b�V���_�����ӏ��j
            YM(I+1) = ZWR(IWR); %! �k���O��LCMS��z���W�iCCS���b�V���_�����ӏ��j
            break
        end
    end
end
XM(1) = RWR(1); %! �k���O��LCMS��r���W�iX�_�ʒu�j
YM(1) = ZWR(1); %! �k���O��LCMS��z���W�iX�_�ʒu�j
for I = 1:NMC
    fprintf(Mesh,'%d   %d\r\n',XM(I),YM(I));
end
figure('Name','Mesh point for plasma current density reconstruction','NumberTitle','off')
plot(XM,YM,'o', BOUND(:,1),BOUND(:,2));
axis equal
%
NNC = 3*NEC;
XM(NMC+1) = XM(1);
YM(NMC+1) = YM(1);
fclose(Mesh);
%-----------------------------------------------------------------------
% �ߓ_���W�̍쐬
XN = zeros(1,3*NEC);
YN = zeros(1,3*NEC);
for I = 1:NEC
    XN(3*I-2) = (5.D0*XM(2*I-1)+5.D0*XM(2*I)-XM(2*I+1))/9.D0;
    YN(3*I-2) = (5.D0*YM(2*I-1)+5.D0*YM(2*I)-YM(2*I+1))/9.D0;
    XN(3*I-1) = XM(2*I);
    YN(3*I-1) = YM(2*I);
    XN(3*I)   = (5.D0*XM(2*I+1)+5.D0*XM(2*I)-XM(2*I-1))/9.D0;
    YN(3*I)   = (5.D0*YM(2*I+1)+5.D0*YM(2*I)-YM(2*I-1))/9.D0;
end
%-----------------------------------------------------------------------
fprintf('%s   %d\n','Number of Elements   = ', NEC);
fprintf('%s   %d\n','Number of Mesh Point = ', NMC);
fprintf('%s   %d\n','Number of Node Point = ', NNC);
%-----------------------------------------------------------------------
%% �ߓ_��̕����ʂ�^����
FI = zeros(1,MXCCS);
DFI = zeros(1,MXCCS);%  ! BOUNDARY CONDITION FI:��, DFI:d��/dn
DFI(1:sum(NCCN)) = FFOUT(1:sum(NCCN));
FI(1:sum(NCCN)) = FFOUT(1 + sum(NCCN):sum(NCCN) + sum(NCCN));
for II = 1:NNC 
% *******************************************************************
    [PSIS,PSISA,PSISB] = QINTER(XN(II),YN(II),GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
    n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,...
    RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS); % OK
% ********************************************************************CC 
    BP = sqrt(PSISA^2 + PSISB^2); 
    fprintf(Node_BC,'%d  %d  %d  %d\r\n',XN(II),YN(II), PSIS, BP);
end
fclose(Node_BC);
end
%% ERROR BOUNDARY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function ERROR_BOUNDARY(IUTST,GETA,RCCS,ZCCS,FFOUT,n100,n50,RC,ZC,ECIGRP,...
ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,...
NONC,fid99,fid100,RMYU0,NE,CCS,MXCCS)
Lout = 0.08;
Lin = 0.08;
REF = dlmread(strcat('UTST_',num2str(IUTST),'/PS.txt'));
% REF(:,1) = spline(1:numel(REFS(:,1)),REFS(:,1),1:1800);
% REF(:,2) = spline(1:numel(REFS(:,2)),REFS(:,2),1:1800);
%REC = dlmread('BOUNDARY_P9n3.TXT');
REC = dlmread(strcat('UTST_',num2str(IUTST),'/PS.txt'));
% REC(:,1) = spline(1:numel(RECS(:,1)),RECS(:,1),1:290);
% REC(:,2) = spline(1:numel(RECS(:,2)),RECS(:,2),1:290);
PP = zeros(numel(REF(:,1)),2);
ZZ = zeros(numel(REF(:,1)),1);
%% �@���x�N�g�������肷��
NVEC = normal_numeric(REF);
REC(end+1,1) = REC(1,1);
REC(end+1,2) = REC(1,2);
for I = 1:numel(REF(:,1))
    for J = 1:numel(REC(:,1))-1
        a1 = REF(I,:);
        if (inpolygon(REF(I,1),REF(I,2), REC(:,1), REC(:,2)) == 1)
            a2 = REF(I,:) - NVEC(I,:)*Lout;
            PorM =  1; 
        else
            a2 = REF(I,:) + NVEC(I,:)*Lin;
            PorM = -1; 
        end
        b1 = REC(J,:);
        b2 = REC(J+1,:);
        a1(3) = 0;
        a2(3) = 0;
        b1(3) = 0;
        b2(3) = 0;
       %%a1,a2��b1,b2��[�_�Ƃ�������̌�������
        if and((sum(cross(a2-a1, b1-a1)) * sum(cross(a2-a1, b2-a1)) < 0 ), ...
               (sum(cross(b2-b1, a1-b1)) * sum(cross(b2-b1, a2-b1)) < 0 ));
         % a1,a2��[�_�Ƃ��������b1,b2��[�_�Ƃ�������̌�_�v�Z
            b = b2-b1;
            d1 = sum(abs(cross(b, a1 - b1)));
            d2 = sum(abs(cross(b, a2 - b1)));
            t = d1 / (d1 + d2);
            P = a1 + (a2-a1) * t;
            PP(I,1) = P(1);
            PP(I,2) = P(2);
            ZZ(I) = PorM*sqrt((a1(1)-P(1))^2+(a1(2)-P(2))^2);
        end
    end
end
[row1,col1] = find(PP(:,1));
PP = PP(row1,[1 2]);
REF = REF(row1,[1 2]);
ZZ = ZZ(row1);
DIST = zeros(numel(PP(:,1)),1);
BP_rec   = zeros(numel(PP(:,1)),1);
BP_ref   = zeros(numel(PP(:,1)),1);
DIST(1) = 0;
for I = 2:numel(PP(:,1))
    DIST(I) = DIST(I-1) + sqrt((REF(I-1,1)-REF(I,1))^2+(REF(I-1,2)-REF(I,2))^2);
end
%
FI = zeros(1,MXCCS);
DFI = zeros(1,MXCCS);%  ! BOUNDARY CONDITION FI:��, DFI:d��/dn
DFI(1:sum(NCCN)) = FFOUT(1:sum(NCCN));
FI(1:sum(NCCN)) = FFOUT(1 + sum(NCCN):sum(NCCN) + sum(NCCN));
for II = 1:numel(PP(:,1))
% *******************************************************************
    [PSIS,PSISA,PSISB] = QINTER(PP(II,1),PP(II,2),GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
    n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,...
    RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS); % OK
% ********************************************************************CC 
    BP_rec(II) = sqrt((PSISA/PP(II,1))^2 + (-PSISB/PP(II,1))^2);
%    BP_rec(II) = sqrt((PSISA/PP(II,2))^2 + (-PSISB/PP(II,1))^2);
end
warning('off','all')
warning
%FLUXref = dlmread(strcat('@UTST_FluxProfile_#',num2str(IUTST),'_m2.txt'));
dir = (['./']);
env = load([dir 'equiv.txt']);
Zref = linspace(env(3),env(4),env(1)); % z_min����z_max�̊Ԃ���env(1)�i-0.9985-0.9985��2033�����j�@
Rref = linspace(env(5),env(6),env(2)); % r_min����r_max�̊Ԃ���env(2)�i-0.10815-0.69400��602�����j
%for II = 1:numel(PP(:,1))
II = 1:numel(PP(:,1));
% % *******************************************************************
%     [PSIS,PSISA,PSISB] = QINTER(REF(II,1),REF(II,2),GETA,RCCS,ZCCS,FFOUT,FI,DFI,...
%     n100,n50,RC,ZC,ECIGRP,ECI,KCMX,RCCN,ZCCN,NCCN,KNE,KNN,REV,ZEV,KSE,KSN,...
%     RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,RMYU0,NE,CCS);% OK
% % ********************************************************************CC 
dr = 0.001;
dz = 0.001;
flux2z = interp2(Rref(1:end),Zref(1:end),FLUXref(1:end,1:end), REF(II,1)+dr, REF(II,2),'spline');
flux1z = interp2(Rref(1:end),Zref(1:end),FLUXref(1:end,1:end), REF(II,1)-dr, REF(II,2),'spline');
flux2r = interp2(Rref(1:end),Zref(1:end),FLUXref(1:end,1:end), REF(II,1), REF(II,2)+dz,'spline');
flux1r = interp2(Rref(1:end),Zref(1:end),FLUXref(1:end,1:end), REF(II,1), REF(II,2)-dz,'spline');
Bz(II)=(flux2r-flux1r)/(2*dr)./REF(II,1)/2/pi;
Br(II)=-(flux2z-flux1z)/(2*dz)./REF(II,1)/2/pi;
BP_ref(II) = sqrt(Br(II).^2 + Bz(II).^2);
%end
dlmwrite(strcat('UTST_',num2str(IUTST),'/errordistance.txt'),horzcat(DIST(1:numel(PP(:,1))),ZZ(1:numel(PP(:,1)))));
figure('Name','Error of boundary','NumberTitle','off')
plot(DIST(1:numel(PP(:,1))),ZZ(1:numel(PP(:,1))),'-r')
xlabel({'Distance [m]'});
ylabel({'Relative error [m]'});
%axis([0 DIST(end) -3*10^(-3) 3*10^(-3)])
set(gca, 'FontSize',14);
dlmwrite(strcat('UTST_',num2str(IUTST),'/errorbpref.txt'),horzcat(DIST(1:numel(PP(:,1))),BP_ref(1:numel(PP(:,1)))));
dlmwrite(strcat('UTST_',num2str(IUTST),'/errorbp.txt'),horzcat(DIST(1:numel(PP(:,1))),BP_rec(1:numel(PP(:,1)))));
figure('Name','Error Bp on boundary','NumberTitle','off')
plot(DIST(1:numel(PP(:,1))),BP_ref(1:numel(PP(:,1))),'--k', DIST(1:numel(PP(:,1))),BP_rec(1:numel(PP(:,1))),'-r')
xlabel({'Distance [m]'});
ylabel({'poloidal magnetic field [T]'});
%axis([0 DIST(end) -3*10^(-3) 3*10^(-3)])
set(gca, 'FontSize',14);
end