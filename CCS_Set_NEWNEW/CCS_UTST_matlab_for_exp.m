function CCS_UTST_matlab()
%% *****************************************************************
%% *  CAUCHY CONDITION SURFACE METHOD  (quadratic element version) *
%% *****************************************************************

format long e;

mu0 =   4.0*pi*1.0e-7;
%MXCCS = 20;          % MAX! / NUMBER OF ELEMENTS ON THE CCS
%MXINT = 10000;       % MAX! / NUMBER OF INTERNAL POINTS


PARAM = loadinputfile;
PARAM.input_file_directory='UTST_exp';

%REF =   loadreference(PARAM);

[SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = loadsensordata(PARAM);

WALL =  loadwalldata(PARAM);

ExtCOIL = loadcoildata(PARAM);
% 
% FFDAT = makeFFdata(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP);
% 
GHR = zeros(1,300); 
GHZ = zeros(1,300); 
CCSDAT = makeCCSdata(PARAM, GHR, GHZ);
% 
% 
% FF = FFDAT;
% FF(end+1:end+sum(CCSDAT.NCCN)) = 0.0;
% if PARAM.IPCONST
%     FF(end+1) = -66410*mu0;
% end
% 
% 
% %% *************************************************************
% %% *************************************************************
% %%      iteration Start
%  
% %CC ITMX=1;
% ITMX = 1000;
% 
% DELGE = 0.0;
% OLDEL = 1.0;
% GETA = 0.0;     %! ���ʂ̏����l���[���ɂ���
% OLDGT = GETA;
% 
% fid90 = fopen([PARAM.temporary_file_directory '\GETA.txt'],'w'); %90
% OMGA = 1.9;
% 
% if (PARAM.ITSKP > 0)
%     ITMX=1;
%     fprintf('���ʃT�[�`�� (1) SVD_MT�����݂̂ł��܂�\n');
% else
%     fprintf('���ʃT�[�`�� (0) INTER�̔����ł����܂�\n');
% end
% 
% ITEND = 0;
% EPS1 = 0;
%    
% for IT = 1:ITMX
%     GETA = GETA + DELGE;
%     GETA = OMGA*GETA + (1.0D0-OMGA)*OLDGT;
% 
%     fprintf('****************************************************\n');
%     fprintf('%s%d%s\n','** (IT) Iteration count =',IT,' **');
%     fprintf('%s%d%d\n','** (GETA/DELGE) =',GETA,DELGE);
% 
%     %******************************************************************************************************
%        
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     %[FC,BR,BZ,PSIFLX,PSIC,AMYU0,AA,FF,GETA,NCCS,fid99,fid100] = FORM(WAHAHA,IT,OLDGT,...
%     %GETA,FF,AA,n100,n50,ECI,KCMX,RS,ZS,RC,ZC,ITYPE,TET,NTPB,NNPB,NAPB,NFLX,...
%     %RCCN,ZCCN,NCCN,NCCS,RCCS,ZCCS,REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,Nedp,NONC,MXCCS,...
%     %RMYU0,NE,CCS,ipconst);%NE uhiski OK
% 
%     AA=zeros(length(FF),150);
% %    AA=zeros(length(FF),sum(CCSDAT.NCCN)*2+WALL.KNN);
%     
%      [FC,BR,BZ,PSIFLX,PSIC,AA,FF,GETA] =...
%      FORM(PARAM,AA,FF,IT,OLDGT,GETA,ExtCOIL,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL);
%  
%     %% INOMOTO start
%     fluxfactor=100;
% 
%     FF(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) =    FF(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) * fluxfactor;
%     AA(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM,:) =  AA(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM,:) * fluxfactor;
%     FC(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) =    FC(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) * fluxfactor;
% 
%     %SENSOR_FLXLP.FLXLP = SENSOR_FLXLP.FLXLP /2/pi * fluxfactor;
%     %FLXLP = FF(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM);
%     FLXLP = FFDAT(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) * fluxfactor;
% 
%     %% INOMOTO end
% 
%     %% Solve Matrix Equation
%     %% Singular Value Decompoition
%     
%     if (PARAM.IPCONST == 1)
%         NMAX = SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN)+1;% ushikiip
%     else
%         NMAX = SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN);% ushikiip
%     end
%     
%     JMAX = sum(CCSDAT.NCCN) + sum(CCSDAT.NCCN) + sum(WALL.KNN) + sum(WALL.KSN);
%     fprintf('NCCN/KNN/KSN = %d %d %d\n',sum(CCSDAT.NCCN),sum(WALL.KNN),sum(WALL.KSN));
%     
%     %  Modified Truncated SVD for Smoothness
%     %[C,W,U,V,FFOUT,XBFR,XMT,XGETA,GET] = SVD_MT(ITSKP,IT,AA,FF,NMAX,JMAX,n100,n50,...
%     %    0,0.0D0,NAPB,NFLX,NCCN,KNN,KSN,FC,FLXLP,BSNSR,GETA_YN,AUTO,AUTOINPUT);
%     [C,W,U,V,FFOUT,XBFR,XMT,XGETA,GET] = SVD_MT_matlab(PARAM,PARAM.ITSKP,IT,AA,FF,FC,NMAX,JMAX,250,150,...
%         0,0.0D0,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL,FLXLP);
%     
%     %���m����
%     half_norm = sqrt((sum((FFOUT).^2)));
%     fprintf('%s%d\r\n','norm of the solution vector = ',half_norm);
% 
%     %% ***************************
%     %% ***************************
%     if (PARAM.ITSKP > 0)
%         GETA = XGETA;
%         break
%     else
%     end
%     
%     %% C***************************
%     %% C***************************
%     %[PSI,DELGE,RCCS,ZCCS,XPSI] = INTER(0,GETA,RCCS,ZCCS,CR,CZ,FFOUT,n100,n50,RC,ZC,...
%     %ECIGRP,ECI,KCMX,NAPB,NFLX,FLXLP,BSNSR,RS,ZS,TET,RCCN,ZCCN,...
%     %NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,WAHAHA,...
%     %MXCCS,MXINT,NCCS,NINT,RMYU0,NE,CCS);
%     %% C***************************
%     %% C***************************
% 
%     EPSDEL = abs((DELGE-OLDEL)/OLDEL);
%     if and(IT > 2, EPSDEL < 0.1)
%         OMGA = OMGA + 5.0;
%         fprintf('%s %d\r\n','    ******  OMGA changed to', OMGA);
%     end
%     
%     OLDEL = DELGE;
%     OLDGT = GETA;
%     EPS1 = abs(DELGE/GETA);
%     ITEND = IT;
%     
%     if (EPS1 < 1.0D-5)
%         break
%     end
%     
% end
% 
% %  ******************************************************************************************************
% %      iteration End
% %  ******************************************************************************************************
% 
% fprintf('********* Congraturation!!  GETA iteration Converged.\n');
% fprintf('%s%d\n','****** No. of iterations required = ',ITEND);
% fprintf('%s%d %d\n','** (GETA/DELGE) = ',GETA,DELGE);
% fprintf('%s%d\n','** EPS=DABS(DELGE/GETA)= ',EPS1);

%% *************************************************************
%% *************************************************************

if 1
    % plot sensor position
    VV = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);
    SEN0 = dlmread([PARAM.temporary_file_directory '/SENPOS0.txt']);
    SEN1 = dlmread([PARAM.temporary_file_directory '/SENPOS1.txt']);
    
    figure('Name','Sensor Position','NumberTitle','off')
    subplot(1,2,1);
    plot(VV(:,1),VV(:,2),'-k');
    hold on
    plot(SEN0(:,1),SEN0(:,2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',8)
    plot(SEN1(:,1),SEN1(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8)
    title('Sensor Position')
    xlabel('r [m]');
    ylabel('z [m]');
    axis equal
    axis([0 1 -1.2 1.2])
    set(gca, 'FontSize',14);

    % plot CV segment
    %COIL = dlmread('output/@UTST_CoilGeom.txt');
    VVMESH = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);
    VVNODE = dlmread([PARAM.temporary_file_directory '/VacuumVesselNodePoints.txt']);
    VVSEG = dlmread([PARAM.temporary_file_directory '/VacuumVesselSegments.txt']);
    subplot(1,2,2);
    %figure('Name','CVsegment','NumberTitle','off')
    plot(VV(:,1),VV(:,2),'-k')
    hold on
    plot(VVMESH(:,1),VVMESH(:,2),'bx','MarkerSize',8)
    hold on
    plot(VVNODE(:,1),VVNODE(:,2),'ko','MarkerFaceColor','k','MarkerSize',4)
    % plot(COIL(:,1),COIL(:,2),'s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
    % plot(RES(1:3),ZES(1:3),'ko',RES(1:3),ZES(1:3),'-k')
    title('CVsegment')
    xlabel('r [m]');
    ylabel('z [m]');
    axis equal
    axis([0 1 -1.2 1.2])
    set(gca, 'FontSize',14);

    for i = 1:PARAM.CCS
        RCCSSP(i,:) = spline(1:CCSDAT.NCCN+1,horzcat(CCSDAT.RCCN(i,1:CCSDAT.NCCN),CCSDAT.RCCN(i,1)),1:1/10:CCSDAT.NCCN+1);
        ZCCSSP(i,:) = spline(1:CCSDAT.NCCN+1,horzcat(CCSDAT.ZCCN(i,1:CCSDAT.NCCN),CCSDAT.ZCCN(i,1)),1:1/10:CCSDAT.NCCN+1);
        
        plot(RCCSSP(i,:),ZCCSSP(i,:),'-k')
        plot(CCSDAT.RCCS(i,1:CCSDAT.NCCS),CCSDAT.ZCCS(i,1:CCSDAT.NCCS),'bx','MarkerSize',8)
        plot(CCSDAT.RCCN(i,1:CCSDAT.NCCN),CCSDAT.ZCCN(i,1:CCSDAT.NCCN),'ko','MarkerFaceColor','k','MarkerSize',4)
    end
end




%% Data positions for 2D plot
PLOT.R = 0.1:0.01:0.9;
PLOT.Z = -1:0.02:1;


end


%% *************************************************************
%% *************************************************************
function PARAM = loadinputfile()

    inputdata0=fileread('CCS_input/input.txt');
    inputdata=strsplit(inputdata0,{'\n','\t','\r'});

    % file direcory
    PARAMNUM=1;
    PARAM.input_file_directory = inputdata{PARAMNUM}; PARAMNUM=PARAMNUM+2;
    
    PARAM.temporary_file_directory = inputdata{3}; PARAMNUM=PARAMNUM+2;
    
    PARAM.output_file_directory = inputdata{5}; PARAMNUM=PARAMNUM+2;

    % CCS parameters
    PARAM.IUTST =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.GETA_YN = str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.NONC =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.SIGM =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.SEED =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.CCS =     str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.IDECCS =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;

    for i=1:PARAM.CCS
        PARAM.R0(i) =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.Z0(i)=    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.RR(i) =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.CAPPER(i) =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.NE(i) =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.MSEC(i,1) =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.MSEC(i,2) =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.MSEC(i,3) =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.SOU(i) =      str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    end
    
    PARAM.IPCONST =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.LCOND =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.ITSKP =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.ITRNC =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.IDCN =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.KUP0 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.MTS =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.LXMAP =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.PLOTW =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.SERX =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.LOCAT =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.AAA1 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.BBB1 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.AAA2 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.BBB2 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.statR =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.endR =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.I0 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.TRUTOT =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
end
%% *************************************************************
%% *************************************************************
function REF = loadreference(PARAM)

    REF.Flux = dlmread([PARAM.input_file_directory '/FluxProfile_2D.txt']);
    REF.Flux = REF.Flux / 2./pi;
    [znum, rnum] = size(REF.Flux);

    %% Here, the geometory should be derived from input file
    zmin = -0.9985;
    zmax = 0.9985;
    rmin = 0.108150;
    rmax = 0.6940;
    delz = (zmax-zmin)/(znum-1);
    delr = (rmax-rmin)/(rnum-1);

    for i=1:znum
        REF.Z(i) = zmin + delz * (i-1);
    end
    for i=1:rnum
        REF.R(i) = rmin + delr * (i-1);
    end
    REF.Flux(abs(REF.Z)>0.285, REF.R > 0.5985) = 0;
end
%% *************************************************************
%% *************************************************************
function [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP] = loadsensordata(PARAM)

    sensordata_B0 = fileread([PARAM.input_file_directory '/Sensor_B.txt']);
    sensordata_B =  strsplit(sensordata_B0,{'\n','\t','\r'});
    sensornum_B =   (length(sensordata_B)-1)/5 - 1;

    chnum = 0;
    for i=1:sensornum_B
        R = str2double(sensordata_B{6+(i-1)*5});
        Z = str2double(sensordata_B{7+(i-1)*5});
        %PSI = str2double(sensordata_B{8+(i-1)*5});
        BZ = str2double(sensordata_B{9+(i-1)*5});
        BR = str2double(sensordata_B{10+(i-1)*5});
        
        if 0%((R< 0.57+0.0015 && R>0.27-0.0015)||(Z<0.239 && R>0.67)||(R<0.689 && R>0.629-0.0015))
        else
            chnum = chnum+1;
            
            %SENSOR_TPRB.R(chnum*2-1) =       R;
            SENSOR_TPRB.R(chnum) =         R;
            %SENSOR_TPRB.Z(chnum*2-1) =       Z;
            SENSOR_TPRB.Z(chnum) =         Z;
            %SENSOR_TPRB.TET(chnum*2-1) =     atan2(BZ,BR);
            SENSOR_TPRB.TET(chnum) =       atan2(BZ,BR);
            %SENSOR_TPRB.TPRB(chnum*2-1) =    sqrt(BR^2 + BZ^2);
            SENSOR_TPRB.TPRB(chnum) =      sqrt(BR^2 + BZ^2);
            %SENSOR_TPRB.ITYPE(chnum*2-1) =   1;
            SENSOR_TPRB.ITYPE(chnum) =     1;
        end
    end
    SENSOR_TPRB.NUM = chnum;
    disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);
 
    %% No NPRB
    SENSOR_NPRB.NUM = 0;
    SENSOR_NPRB.R = [];
    SENSOR_NPRB.Z = [];
    SENSOR_NPRB.TET = [];
    SENSOR_NPRB.NPRB = [];
    SENSOR_NPRB.ITYPE = [];
    
    sensordata_Flux0 =  fileread([PARAM.input_file_directory '/Sensor_Flux.txt']);
    sensordata_Flux =   strsplit(sensordata_Flux0,{'\n','\t','\r'});
    sensornum_Flux =    (length(sensordata_Flux)-1)/5 - 1;
    
    chnum = 0;
    
    for i=1:sensornum_Flux
        R = str2double(sensordata_Flux{6+(i-1)*5});
        Z = str2double(sensordata_Flux{7+(i-1)*5});
        PSI = str2double(sensordata_Flux{8+(i-1)*5});
        %BZ = str2double(sensordata_Flux{8+(i-1)*5});
        %BR = str2double(sensordata_Flux{10+(i-1)*5});

        chnum = chnum+1;

        %SENSOR_FLXLP.R(chnum*2-1) =       R;
        SENSOR_FLXLP.R(chnum) =         R;
        %SENSOR_FLXLP.Z(chnum*2-1) =       Z;
        SENSOR_FLXLP.Z(chnum) =         Z;
        %SENSOR_FLXLP.FLXLP(chnum*2-1) =   PSI;
        SENSOR_FLXLP.FLXLP(chnum*2) =     PSI;
        %SENSOR_FLXLP.TET(chnum*2-1) =     0.0D0;
	    SENSOR_FLXLP.TET(chnum*2) =       0.0D0;
	    %SENSOR_FLXLP.ITYPE(chnum*2-1) =   0;
	    SENSOR_FLXLP.ITYPE(chnum*2) =     0;
    end
    
    SENSOR_FLXLP.NUM=chnum;
    disp(['Number of FLXLP = ' num2str(SENSOR_FLXLP.NUM)]);

    %% Write files for CCS
    fp = fopen([PARAM.temporary_file_directory '/SENPOS0.txt'],'w'); % 110
    for i=1:SENSOR_FLXLP.NUM
        fprintf(fp,'%d %d\n',SENSOR_FLXLP.R(i),SENSOR_FLXLP.Z(i));
    end
    fclose(fp);
    
    fp = fopen([PARAM.temporary_file_directory '/SENPOS1.txt'],'w'); % 111
    for i=1:SENSOR_TPRB.NUM
        fprintf(fp,'%d %d\n',SENSOR_TPRB.R(i),SENSOR_TPRB.Z(i));
    end
    fclose(fp);
    
    fp = fopen([PARAM.temporary_file_directory '/SENPOS2.txt'],'w'); % 112
    for i=1:SENSOR_NPRB.NUM
        fprintf(fp,'%d %d\n',SENSOR_NPRB.R(i),SENSOR_NPRB.Z(i));
    end
    fclose(fp);
 
    %%  *************************************************************************
    %%     Generation of CCS input data
    %%  *************************************************************************

    fprintf('Generation of CCS input data ** START ***\n');
    fp = fopen([PARAM.temporary_file_directory '/CCSinput_UTST(temp).txt'],'w');
    
    fprintf(fp,'%s\n','*');
    fprintf(fp,'%s\n','*** CCS Test input for UTST generated in PreUTST ***'); 
    fprintf(fp,'%s\n','** NTPB/NNPB/NFLX=(No. of T-/N-probes/Flux-Loops) **');
    fprintf(fp,'   %d     %d     %d\n',SENSOR_TPRB.NUM,SENSOR_NPRB.NUM,SENSOR_FLXLP.NUM); 
    fprintf(fp,'%s\n','**  GETA (SSURF)');
    fprintf(fp,'  %s\n','0.0000E+00');
    fprintf(fp,'%s\n','* T-Probe');
    for i = 1:SENSOR_TPRB.NUM
        fprintf(fp,' %d\n',SENSOR_TPRB.TPRB(i));
    end
    fprintf(fp,'%s\n','* N-Probe');
    fprintf(fp,'%s\n','* Flux-Loop');
    for i = 1:SENSOR_FLXLP.NUM
        fprintf(fp,' %d\n',SENSOR_FLXLP.FLXLP(i));
    end
    fprintf(fp,'%s\n','****** MINR * MAXR * MINZ * MAXZ ****');
    fprintf(fp,'%s\n','10   90  -100  100');
    fprintf(fp,'%s\n','*********');
    fprintf(fp,'%s\n','* ---�R�C���d���f�[�^�̕���---�P��[kA]');
    fprintf(fp,'%s\n','* EF');
    fprintf(fp,'%s\n','* PF#1');
    fprintf(fp,'%s\n','* PF#2');
    fprintf(fp,'%s\n','* PF#3');
    fprintf(fp,'%s\n','* PF#4');

    fprintf(fp,'   %s\n','-0.28');
    fprintf(fp,'   %s\n','0.0');
    fprintf(fp,'   %s\n','0.0');
    fprintf(fp,'   %s\n','0.0');
    fprintf(fp,'   %s\n','0.0');
    
    fclose(fp);
    fprintf('Generation of CCS input data ** END ***\n');
end
%% *************************************************************
%% *************************************************************
function WALL =  loadwalldata(PARAM)

    %% Here, the geometory should be derived from input file
%     R = [0.694D0, 0.694D0, 0.5985D0, 0.5985D0, 0.10815D0, 0.10815D0,... 
%         0.10815D0, 0.5985D0, 0.5985D0, 0.694D0,  0.694D0];
%     Z = [0.0D0, 0.285D0, 0.285D0, 0.9985D0, 0.9985D0, 0.0D0,...
%         -0.9985D0, -0.9985D0, -0.285D0, -0.285D0 , 0.0D0];     


%% CONDSHL_UTST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ���̏�̉Q�d���ߓ_�ʒu�f�[�^�̐ݒ�Ɠǂݎ��      
%function [KNE,KNN,REV,ZEV,RES,ZES,KSE,KSN,REVN,ZEVN] = CONDSHL_UTST(WAHAHA,Nedp,NONC)

    MAXM = 13;
    Z = 0.15;
    INBOAD = 0.5;

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

    fp = fopen([PARAM.temporary_file_directory '/VacuumVesselSegments.txt'],'w');
    for i=1:MAXM+1
        fprintf(fp,'%d %d\n',RSEC(i),ZSEC(i));
    end
    fclose(fp);

    % �^��e���̉Q�d���ߓ_�̐ݒ�
    % KNE=�^��e��̕������E�v�f��,  KNM=�^��e���̃��b�V���_��,  KNN=�^��e���̐ߓ_��

    II=1;
    WALL.REV(II)=RSEC(1);
    WALL.ZEV(II)=ZSEC(1);

    WALL.KNE=0;
    for I = 1:MAXM
        WALL.KNE = WALL.KNE+MSEC(I);
        CM = 2*MSEC(I);
        DELR = (RSEC(I+1)-RSEC(I))/CM;
        DELZ = (ZSEC(I+1)-ZSEC(I))/CM;
        for J = 1:2*MSEC(I)
            II = II+1;
            WALL.REV(II) = WALL.REV(II-1)+DELR;
            WALL.ZEV(II) = WALL.ZEV(II-1)+DELZ;
        end
    end

    KNM = II-1;

    fp = fopen([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt'],'w');
    for II = 1:KNM+1
        fprintf(fp,'%d %d\n',WALL.REV(II),WALL.ZEV(II));
    end
    fclose(fp);

    %% �I�I�I�I�I�I�I�@ ��K���v�f�ߓ_�v�Z�̗\��n�@�@�@�I�I�I�I�I�I�I
    % �^��e���̔�K���v�f�ߓ_���W�̍쐬
    if (PARAM.NONC > 0)
        WALL.KNN = WALL.KNE*3;
        I = 1:WALL.KNE;
	    WALL.REVN(3*I-2) = (5.D0*WALL.REV(2*I-1)+5.D0.*WALL.REV(2*I)-WALL.REV(2*I+1))/9.D0;
	    WALL.ZEVN(3*I-2) = (5.D0*WALL.ZEV(2*I-1)+5.D0.*WALL.ZEV(2*I)-WALL.ZEV(2*I+1))/9.D0;
	    WALL.REVN(3*I-1) = WALL.REV(2*I);
	    WALL.ZEVN(3*I-1) = WALL.ZEV(2*I);
	    WALL.REVN(3*I)   = (5.D0*WALL.REV(2*I+1)+5.D0.*WALL.REV(2*I)-WALL.REV(2*I-1))/9.D0;
	    WALL.ZEVN(3*I)   = (5.D0*WALL.ZEV(2*I+1)+5.D0.*WALL.ZEV(2*I)-WALL.ZEV(2*I-1))/9.D0;
    else
        WALL.KNN = KNM;
    end

    fp = fopen([PARAM.temporary_file_directory '/VacuumVesselNodePoints.txt'],'w'); %113
    if (PARAM.NONC == 0)
        for II = 1:KNM+1
            fprintf(fp,'%d %d\n',WALL.REV(II),WALL.ZEV(II));
        end
        fprintf('KNE=%d    KNM=%d\n',WALL.KNE,KNM);
        fprintf('%s\n','The end meshpoint agrees with the start meshpoint');
    else
        for II = 1:WALL.KNN
            fprintf(fp,'%d %d\n',WALL.REVN(II),WALL.ZEVN(II));
        end
        fprintf('KNE/KNM/KNN = %d  %d  %d\n',WALL.KNE,KNM,WALL.KNN);
    end
    fclose(fp);

    % ���艻��̉Q�d���ߓ_�̐ݒ�
    % KSE=���艻�̕������E�v�f��,  KSN=���艻��̐ߓ_��
    %***
    %cd      KSE=1 
    %prompt = '�d���V�[�g��ɋ��E�v�f��z�u����H (Yes/No)=(1/0)\n';
    %KSE = input(prompt);
    
    WALL.KSE = 0;
    %***
    WALL.KSN = WALL.KSE*2 + 1;
    if (WALL.KSE == 0)
        WALL.KSN=0;
    end
    if (WALL.KSE > 0)
        fp = fopen([PARAM.temporary_file_directory /StabilizerPoints.txt'],'w'); % 114
        WALL.RES(1) = 0.25;
        WALL.RES(2) = 0.28;
        WALL.RES(3) = 0.31;
        WALL.ZES(1) = 0;     
        WALL.ZES(2) = 0;    
        WALL.ZES(3) = 0;     
        for I = 1:WALL.KSN
            fprintf(fp,'%d %d\n',WALL.RES(I),WALL.ZES(I));
        end
        fclose(fp);
    else
        WALL.RES = [];
        WALL.ZES = [];     
    end
end
%% *************************************************************
%% *************************************************************
function ExtCOIL = loadcoildata(PARAM)

    %% Here, the geometory should be derived from input file
    ExtCOIL.NUM = 10;
    ExtCOIL.NAME = ['EFL', 'EFU', 'PF1L', 'PF1U', 'PF2L', 'PF2U',  'PF3L', 'PF3U', 'PF4L', 'PF4U'];
    ExtCOIL.R = [0.80, 0.80, 0.20, 0.20, 0.665, 0.665, 0.750, 0.750, 0.685, 0.685];
    ExtCOIL.Z = [-1.07, 1.07, -1.10, 1.10, -0.80, 0.80, -0.675, 0.675, -0.50, 0.50];
    ExtCOIL.C = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    ExtCOIL.N = [200, 200, 8, 8, 3, 3, 8, 8, 3, 3]; 
    ExtCOIL.I = [0.28, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
end
%% *************************************************************
%% *************************************************************
function FFDAT = makeFFdata(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP)

    %% FFDAT0 (w/o noise)
    for i=1:SENSOR_TPRB.NUM
        FFDAT0(i)=SENSOR_TPRB.TPRB(i);
    end
    for i=1:SENSOR_NPRB.NUM
        FFDAT0(i+SENSOR_TPRB.NUM)=SENSOR_NPRB.NPRB(i);
    end
    for i=1:SENSOR_FLXLP.NUM
        FFDAT0(i+SENSOR_TPRB.NUM+SENSOR_NPRB.NUM)=SENSOR_FLXLP.FLXLP(i)/2/pi;
    end

    %% FFDAT (w. noise)
    rng(PARAM.SEED);
    GASDEV = randn(SENSOR_TPRB.NUM+SENSOR_NPRB.NUM+SENSOR_FLXLP.NUM,1);
    for i=1:SENSOR_TPRB.NUM+SENSOR_NPRB.NUM+SENSOR_FLXLP.NUM
        FFDAT(i) = FFDAT0(i) * (1.0 + PARAM.SIGM*GASDEV(i,1));
    end

end
%% *************************************************************
%% *************************************************************
function CCSDAT = makeCCSdata(PARAM, GHR, GHZ)

    PARAM.NE = PARAM.MSEC(1:PARAM.CCS,1) + PARAM.MSEC(1:PARAM.CCS,2) + PARAM.MSEC(1:PARAM.CCS,3);
    %PARAM.NCCS = NE*2;

    ICONT = 0;
    %  9999 CONTINUE ! ICONT�̍X�V�iCCS����蒼���j!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (ICONT == 0)  % ICONT>0 �Ȃ�CCS����蒼���B
        if (PARAM.IDECCS > 0) %! �ȉ~�^CCS or D�^CCS
            %*******************************************************
            %  'D'�^CCS                                         ****     
            %*******************************************************
            CCSDAT = D_CCS(PARAM);
        else
            % *******************************************************
            %   �ȉ~�^CCS                                        ****     
            % *******************************************************
            CCSDAT.NCCS = PARAM.NE*2;
            CCSDAT.NCCN = PARAM.NE*3;

            for i=1:PARAM.CCS
                DTHETA = 2.0*pi/CCSDAT.NCCS(i);
                TRIG = 0.0;      %�O�p�x  
                for j=1:CCSDAT.NCCS(i)
                    THETA = pi/2.0-DTHETA*(j-1);  % CCS�ł̐ϕ��́A���v���
                    CCSDAT.RCCS(i,j) = PARAM.R0(i) + PARAM.RR(i)*cos(THETA + asin(TRIG)*sin(THETA));
                    CCSDAT.ZCCS(i,j) = PARAM.Z0(i) + PARAM.CAPPER(i)*PARAM.RR(i)*sin(THETA);
                end
            end
        end
    else
        % *******************************************************
        %    LCMS�ɑ�����CCS    (ICONT>0 �̂Ƃ�)              ****
        % *******************************************************
        for i=1:PARAM.CCS
            for j=1:CCSDAT.NCCS(i)
                k=CCSDAT.NCCS(i)+2-j;
                if (k>NCCS(i))
                    k=1;
                end
                CCSDAT.RCCS(i,j) = GHR(k);
                CCSDAT.ZCCS(i,j) = GHZ(k);
            end 
        end
    end

    % *******************************************************
    %  CCS��̔�K���v�f�ߓ_���W�̍쐬
    % *******************************************************
    fid12 = fopen([PARAM.temporary_file_directory '/MeshPoints.txt'],'w'); 
    fid13 = fopen([PARAM.temporary_file_directory '/DiscontinuousNodePoints.txt'],'w');
    fid14 = fopen([PARAM.temporary_file_directory '/CurvCCS.txt'],'w');
    fid15 = fopen([PARAM.temporary_file_directory '/CurvCCS_Final.txt'],'w'); 
    
    for i = 1:PARAM.CCS
        CCSDAT.RCCS(i,CCSDAT.NCCS+1) = CCSDAT.RCCS(i,1);
        CCSDAT.ZCCS(i,CCSDAT.NCCS+1) = CCSDAT.ZCCS(i,1);
        
        I = 1:PARAM.NE(i);
        
        CCSDAT.RCCN(i,3*I-2) = (5.*CCSDAT.RCCS(i,2*I-1) + 5.*CCSDAT.RCCS(i,2*I) - CCSDAT.RCCS(i,2*I+1))/9;
        CCSDAT.ZCCN(i,3*I-2) = (5.*CCSDAT.ZCCS(i,2*I-1) + 5.*CCSDAT.ZCCS(i,2*I) - CCSDAT.ZCCS(i,2*I+1))/9;
        CCSDAT.RCCN(i,3*I-1) = CCSDAT.RCCS(i,2*I);
        CCSDAT.ZCCN(i,3*I-1) = CCSDAT.ZCCS(i,2*I);
        CCSDAT.RCCN(i,3*I)   = (5.*CCSDAT.RCCS(i,2*I+1) + 5.*CCSDAT.RCCS(i,2*I) - CCSDAT.RCCS(i,2*I-1))/9;
        CCSDAT.ZCCN(i,3*I)   = (5.*CCSDAT.ZCCS(i,2*I+1) + 5.*CCSDAT.ZCCS(i,2*I) - CCSDAT.ZCCS(i,2*I-1))/9;
        
        for j=1:CCSDAT.NCCS(i)
            fprintf(fid12,'%d %d\n',CCSDAT.RCCS(i,j),CCSDAT.ZCCS(i,j));
        end
        fprintf(fid12,'%d %d\n',CCSDAT.RCCS(i,1),CCSDAT.ZCCS(i,1));
        
        for j=1:CCSDAT.NCCN(i)
            fprintf(fid13,'%d %d\n',CCSDAT.RCCN(i,j),CCSDAT.ZCCN(i,j)); 
        end
        
        %  Draw a fine curve to express the shape of CCS 
        MAXG = 200;
        GMAX = MAXG;
        DEL =   2.0/GMAX;
        for I = 1:PARAM.NE(i)
            for J = 1:MAXG+1
                CJM1 = J-1;
                GII = -1.0+DEL*CJM1;
                F1 = GII*(GII-1.0D0)*0.5;
                F2 = 1.0D0-GII^2;
                F3 = GII*(GII+1.0D0)*0.5;
                RGI = CCSDAT.RCCS(i,2*I-1)*F1 + CCSDAT.RCCS(i,2*I)*F2+CCSDAT.RCCS(i,2*I+1)*F3;
                ZGI = CCSDAT.ZCCS(i,2*I-1)*F1 + CCSDAT.ZCCS(i,2*I)*F2+CCSDAT.ZCCS(i,2*I+1)*F3;
                
                fprintf(fid14,'%d %d\n',RGI,ZGI); 
                fprintf(fid15,'%d %d\n',RGI,ZGI); 
            end
        end
        
    end
    
    fclose(fid12);
    fclose(fid13);
    fclose(fid14);
    fclose(fid15);
    
end
%% *************************************************************
%% *************************************************************
function CCSDAT = D_CCS(PARAM)

    SR(1:PARAM.CCS,1) = PARAM.R0 - PARAM.RR;
    SR(1:PARAM.CCS,3) = SR(1:PARAM.CCS,1);
    SR(1:PARAM.CCS,2) = PARAM.R0 + PARAM.RR;
    SZ(1:PARAM.CCS,1) = PARAM.Z0 + PARAM.RR.*PARAM.CAPPER;
    SZ(1:PARAM.CCS,2) = PARAM.Z0;
    SZ(1:PARAM.CCS,3) = PARAM.Z0 - PARAM.RR.*PARAM.CAPPER;

    fid58 = fopen([PARAM.temporary_file_directory '/@D_CCS_SRpoints.txt'],'w');
    fid59 = fopen([PARAM.temporary_file_directory '/@D_CCS_CheckWrite.txt'],'w');
    fid13 = fopen([PARAM.temporary_file_directory '/DiscontinuousNodePoints.txt'],'w');

    for i=1:PARAM.CCS
        fprintf(fid58,'%d %d\n',PARAM.R0(i),PARAM.Z0(i));
        for j=1:3
            fprintf(fid58,'%d %d\n', SR(i,j),SZ(i,j));
        end
        
        %% �Ȑ����̃��b�V���_���`
        II=0;
        for j=1:2
            M2 = PARAM.MSEC(i,j)*2; % ���b�V���_�̐�
            DEL = 1.0/M2;
            GISTAT = -1.0+(j-1); % �O�U�C -1 0 1
            for k=1:M2
                GI = GISTAT+DEL*(k-1); %�O�U�C �����E�v�f���ɕ������Ă���
                % Compute the values of the shape functions at the integration points
                F1 = GI*(GI-1.0)/2;
                F2 = 1.0-GI^2;
                F3 = GI*(GI+1.0)/2;
                II = II+1;
                CCSDAT.RCCS(i,II) = SR(i,1)*F1+SR(i,2)*F2+SR(i,3)*F3;
                CCSDAT.ZCCS(i,II) = SZ(i,1)*F1+SZ(i,2)*F2+SZ(i,3)*F3; 
            end
        end
        %% �������̃��b�V���_���`
        II=II+1;
        CCSDAT.RCCS(i,II) = SR(i,3);
        CCSDAT.ZCCS(i,II) = SZ(i,3);
        KADO = II;
        M2 = PARAM.MSEC(i,3)*2;
        DELR = (SR(i,1)-SR(i,3))/M2;
        DELZ = (SZ(i,1)-SZ(i,3))/M2;

        for k=1:M2
            II=II+1;
            CCSDAT.RCCS(i,II) = CCSDAT.RCCS(i,II-1) + DELR;
            CCSDAT.ZCCS(i,II) = CCSDAT.ZCCS(i,II-1) + DELZ;
        end
        CCSDAT.NCCS(i) = II-1;
        CCSDAT.NCCN(i) = PARAM.NE(i)*3;

        fprintf('%d%s%d  %d  %d\n',i,'�Ԗ�CCS: NE/NCCS/NCCN = ', PARAM.NE(i), CCSDAT.NCCS(i),CCSDAT.NCCN(i));

        for j=1:CCSDAT.NCCS(i)+1
            fprintf('%d %d %d\n',j, CCSDAT.RCCS(i,j), CCSDAT.ZCCS(i,j));
            fprintf(fid59,'%d %d\n', CCSDAT.RCCS(i,j), CCSDAT.ZCCS(i,j));
            if (II==KADO)
                fprintf(fid59,'%d %d\n',CCSDAT.RCCS(i,j), CCSDAT.ZCCS(i,j));
            end
        end
        %% CCS��̔�K���v�f�ߓ_���W�̍쐬
        CCSDAT.RCCS(i,CCSDAT.NCCS+1) = CCSDAT.RCCS(i,1);
        CCSDAT.ZCCS(i,CCSDAT.NCCS+1) = CCSDAT.ZCCS(i,1);
        CCSDAT.NCCN(i) = PARAM.NE(i)*3; % ��K���v�f
        
        j=1:PARAM.NE(i);
        CCSDAT.RCCN(i,3*j-2) = (5*CCSDAT.RCCS(i,2*j-1) + 5.*CCSDAT.RCCS(i,2*j) - CCSDAT.RCCS(i,2*j+1))/9;
        CCSDAT.ZCCN(i,3*j-2) = (5*CCSDAT.ZCCS(i,2*j-1) + 5.*CCSDAT.ZCCS(i,2*j) - CCSDAT.ZCCS(i,2*j+1))/9;
        CCSDAT.RCCN(i,3*j-1) = CCSDAT.RCCS(i,2*j);
        CCSDAT.ZCCN(i,3*j-1) = CCSDAT.ZCCS(i,2*j);
        CCSDAT.RCCN(i,3*j) = (5*CCSDAT.RCCS(i,2+j+1) + 5.*CCSDAT.RCCS(i,2*j) - CCSDAT.RCCS(i,2*j-1))/9;
        CCSDAT.ZCCN(i,3*j) = (5*CCSDAT.ZCCS(i,2*j+1) + 5.*CCSDAT.ZCCS(i,2*j) - CCSDAT.ZCCS(i,2*j-1))/9;
        
        for j=1:CCSDAT.NCCN(i)
            fprintf(fid13,'%d %d\n',CCSDAT.RCCN(i,j), CCSDAT.ZCCN(i,j));
        end
    end
    fclose(fid13);
    fclose(fid58);
    fclose(fid59);
end


%% FORM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% AA*X=FF
%function [FC,BR,BZ,PSIFLX,PSIC,AMYU0,AA,FF,GETA,NCCS,fid99,fid100]= FORM(WAHAHA,...
%IT,OLDGT,GETA,FF,AA,n100,n50,ECI,KCMX,RS,ZS,RC,ZC,ITYPE,TET,NTPB,NNPB,...
%NAPB,NFLX,RCCN,ZCCN,NCCN,NCCS,RCCS,ZCCS,REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,Nedp,...
%NONC,MXCCS,RMYU0,NE,CCS,ipconst)

function [FC,BR,BZ,PSIFLX,PSIC,AA,FF,GETA] =...
    FORM(PARAM,AA,FF,IT,OLDGT,GETA,ExtCOIL,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL)

    ECI=ExtCOIL.I.*ExtCOIL.N*1000;
    KCMX=ExtCOIL.NUM;
    RS=[SENSOR_FLXLP.R SENSOR_TPRB.R SENSOR_NPRB.R];
    ZS=[SENSOR_FLXLP.Z SENSOR_TPRB.Z SENSOR_NPRB.Z];
    RC=ExtCOIL.R;
    ZC=ExtCOIL.Z;
    %ITYPE=[SENSOR_FLXLP.ITYPE SENSOR_TPRB.ITYPE SENSOR_NPRB.ITYPE];
    TET=[SENSOR_FLXLP.TET SENSOR_TPRB.TET SENSOR_NPRB.TET];
    %NTPB=SENSOR_TPRB.NUM;
    %NNPB=SENSOR_NPRB.NUM;
    NAPB=SENSOR_TPRB.NUM+SENSOR_NPRB.NUM;
    NFLX=SENSOR_FLXLP.NUM;
    RCCN=CCSDAT.RCCN;
    ZCCN=CCSDAT.ZCCN;
    NCCN=CCSDAT.NCCN;
    NCCS=CCSDAT.NCCS;
    RCCS=CCSDAT.RCCS;
    ZCCS=CCSDAT.ZCCS;
    REV=WALL.REV;
    ZEV=WALL.ZEV;
    KNE=WALL.KNE;
    KNN=WALL.KNN;
    RES=WALL.RES;
    ZES=WALL.ZES;
    KSE=WALL.KSE;
    KSN=WALL.KSN;
    %Nedp=150;
    NONC=PARAM.NONC;
    %MXCCS=20;
    RMYU0=4*pi*1e-7;
    NE=PARAM.NE;
    CCS=PARAM.CCS;
    ipconst=PARAM.IPCONST;


     RNOR = 0; 
     ZNOR = 0;
%     BR = zeros(1,NAPB);
%     BZ = zeros(1,NAPB);
%     PSIFLX = zeros(1,NFLX);
%     PSIC = zeros(1,MXCCS);
%     FC = zeros(1,n100);
% 	
%     for I = 1:NAPB+NFLX	
%         fprintf(WAHAHA,'%d Type=%d RS/ZS= %d %d Tet= %d\r\n',I,ITYPE(I),RS(I),ZS(I),TET(I));
%     end

% **********************************************************************
if (IT > 1)
else
% **********************************************************************
%����Ō��̈ʒu�ɖ߂�
% for I = 1:CCS
%     RCCS(I,NCCS(I)+1) = RCCS(I,1);
%     ZCCS(I,NCCS(I)+1) = ZCCS(I,1);
% end
%    !======================================================================
%    !
%    !    FFFFFF
%    !    FF                             | T-PROBE |
%    !    FFFF                VECTOR  FF=| N-PROBE |
%    !    FF                             |FLUX-LOOP|
%    !    FF                             |   CCS   |
%    !
%
%     fprintf(WAHAHA,'%s\r\n','*****�Z���T�[�M������R�C���d����^�̍�����*****');
%  �����̉�FF���쐬�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I�I
%    !FLUX-LOOP 
%     fprintf(WAHAHA,'%s\r\n','FLUX-LOOP     PSI  caused by external coils');
    
    for L = 1:SENSOR_FLXLP.NUM
        %PSIFLX(L) = 0.0D0;
        PSIFLX(L) = GETA;
        
        [PPSIFLX(1:KCMX),PHIR,PHIZ,PPSIA,PPSIB,PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,DAG,DBG,DAH,DBH] =...
            STARB(0,RS(L),ZS(L),RC(1:KCMX),ZC(1:KCMX),RNOR,ZNOR); % OK
	    
        PSIFLX(L) = PSIFLX(L) + sum(PPSIFLX(1:KCMX).*ECI(1:KCMX)*RMYU0);
%        fprintf(WAHAHA,'%d %d\r\n',L,PSIFLX(L));
        FF(L+NAPB) = FF(L+NAPB) - PSIFLX(L);                %! ���ʏ������܂�
        FC(L+NAPB) = PSIFLX(L);                             %! �R�C���d����^
    end
%    !T-PROBE & N-PROBE
%    fprintf(WAHAHA,'%s\r\n','T-PROBE & N-PROBE   B  caused by external coils');

    for L = 1:NAPB
        BR(L) = 0.0D0;
        BZ(L) = 0.0D0;
        
        [PHI,PHIR,PHIZ,PHIA(1:KCMX),PHIB(1:KCMX),PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,DAG,DBG,DAH,DBH] =...
            STARB(1,RS(L+NFLX),ZS(L+NFLX),RC(1:KCMX),ZC(1:KCMX),RNOR,ZNOR); % OK
        
        BR(L) = sum((-PHIB(1:KCMX)/RS(L+NFLX)).*ECI(1:KCMX)*RMYU0);
	    BZ(L) = sum((PHIA(1:KCMX)/RS(L+NFLX)).*ECI(1:KCMX)*RMYU0);
        BBB = BR(L)*cos(TET(L+NFLX))+BZ(L)*sin(TET(L+NFLX));
%        fprintf(WAHAHA,'%d %d\r\n',L,BBB);
        FF(L) = FF(L) - BBB;
        FC(L) = BBB;                                        %! �R�C���d����^           
    end
%    !CCS
%    fprintf(WAHAHA,'%s\r\n','CCS       PSI  caused by external coils');
    for III = 1:CCS
        for L = 1:NCCN(III)
            PSIC(L) = GETA;
            
            [PPSIC(1:KCMX),PHIR,PHIZ,PPSIA,PPSIB,PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,DAG,DBG,DAH,DBH] =...
                STARB(0,RCCN(III,L),ZCCN(III,L),RC(1:KCMX),ZC(1:KCMX),RNOR,ZNOR); % OK

            PSIC(L) = PSIC(L) + sum(PPSIC(1:KCMX).*ECI(1:KCMX)*RMYU0);
            FF(L+sum(NCCN(1:III-1))+NAPB+NFLX) = - PSIC(L); %! ���ʏ������܂ށ@�@
            FC(L+sum(NCCN(1:III-1))+NAPB+NFLX) = PSIC(L);   %! �R�C���d����^�@�@
%            fprintf(WAHAHA,'%d %d %d\r\n',L,PSIC(L),FF(L+sum(NCCN(1:III-1))+NAPB+NFLX));
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
    fid99 = fopen([PARAM.temporary_file_directory '/MINDIST.txt'],'w'); %99
    frewind(fid99);
    fid100 = fopen([PARAM.temporary_file_directory '/SEKIBUNCHECK.PRI'],'w'); %100 
    frewind(fid100);
    fprintf(fid99, '%s\n','****************************************************');
    fprintf(fid99, '%s\n','***    In the Subr. FORM ***************************');
    fprintf(fid99, '%s\n','****************************************************');
    fprintf(fid100,'%s\n','***************************************************');
    fprintf(fid100,'%s\n','***    In the Subr. FORM **************************');
    fprintf(fid100,'%s\n','***************************************************');
%%%%
%!FLUX-LOOP
    for III = 1:CCS
        for I = 1:NFLX
            for K = 1:NE(III)
                [HW,GW,GR,GZ,HR,HZ] = INTEGS(RS(I),ZS(I),RCCS(III,2*K-1),ZCCS(III,2*K-1),...
                    RCCS(III,2*K),ZCCS(III,2*K),RCCS(III,2*K+1),ZCCS(III,2*K+1));% OK
                for JJ = 1:3
                    KK = 3*(K-1)+JJ;
                    AA(I+NAPB,KK+3*sum(NE(1:III-1))) = AA(I+NAPB,KK+3*sum(NE(1:III-1)))+GW(JJ);
                    AA(I+NAPB,KK+3*sum(NE(1:III-1))+sum(NCCN)) = AA(I+NAPB,KK+3*sum(NE(1:III-1))+sum(NCCN))-HW(JJ);
                end
            end% 119
        end
%C      
%!T-PROBE & N-PROBE
        for I = 1:NAPB
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
                if and((3*K) >= I , I >= (3*K-2))
                    %//// ���ْl�ϕ� ///////////////////////////////////////////////////////
                    if (I == (3*K-2)) % ���E��1�Ԗڂ̐ߓ_
                        NODO = 1;
                    else
%                    if(I == (3*K-1))
                        NODO = 2.*(I == (3*K-1)) +3.*(I ~= (3*K-1));% ���E��2�A3�Ԗڂ̐ߓ_
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
            fprintf('%s\n','�^��e���̉Q�d���� flux loop �Z���T�[�ɍ�郵');
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
                        else
                            AA(I+NAPB,sum(NCCN)*2+K-1+JJ) = AA(I+NAPB,sum(NCCN)*2+K-1+JJ)+EE;%  ! ### 2
                        end
                    end
                end
            end
%
            fprintf('%s\n','�^��e���̉Q�d���� ����Z���T�[�ɍ��a');
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
            fprintf('%s\n','�^��e���̉Q�d���� CCS�ɍ�郵');
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
            fprintf('%s\n','�^��e���̉Q�d���� flux loop �Z���T�[�ɍ�郵');
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
                    end
                end
            end
%
            fprintf('%s\n','�^��e���̉Q�d���� ����Z���T�[�ɍ��a');
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
            fprintf('%s\n','�^��e���̉Q�d���� CCS�ɍ�郵');
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
    % ??????????????????????????????????????????????????????????????????    
    %  �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��   �Q�d��
    %    �@���艻�@�@�@�@�@���艻�@�@�@�@�@���艻�@�@�@�@���艻��
    %CAUTION!! The stabilizer is not closed in the poloidal direction.
    %     KSE=No. of boundary elements along the stabilizer
    %     KSN=No. of nodes along the stabilizer (KSN=KNE*2+1)
    %     (RES(),ZES())=Eddy Current Nodes on the stabilizer
    % ??????????????????????????????????????????????????????????????????    
    if (KSE <= 0) % 991 to 992
    else
        AMYU0 = RMYU0*1.0D06;  % ! NAMUAMUdabutsu  #2
        fprintf('%s\n','���艻��̉Q�d���� flux loop �Z���T�[�ɍ�郵');
        for I = 1:NFLX
            A = RS(I);
            B = ZS(I);
            for K = 1:2:KSN-2
                [GW,GR,GZ] =  EXTINDC(A,B,RES(K),ZES(K),RES(K+1),ZES(K+1),RES(K+2),...
                ZES(K+2),NONC,fid99,fid100);% OK 
                for JJ = 1:3
                    EE = GW(JJ)*AMYU0;    %  ! flux*AMYU0         
                    AA(I+NAPB,sum(NCCN)*2+KNN+K-1+JJ) = AA(I+NAPB,sum(NCCN)*2+KNN+K-1+JJ)+EE; % ! ### 2
                end
            end
        end
%CC
        fprintf('%s\n','���艻��̉Q�d���� ����Z���T�[�ɍ��a');
        for I = 1:NAPB
            COST = cos(TET(I+NFLX));
            SINT = sin(TET(I+NFLX));
            A = RS(I+NFLX);
            B = ZS(I+NFLX);
            for K = 1:2:KSN-2
                [GW,GR,GZ] =  EXTINDC(A,B,RES(K),ZES(K),RES(K+1),ZES(K+1),RES(K+2),...
                ZES(K+2),NONC,fid99,fid100); % OK 
                for JJ = 1:3
                    EE = (-COST*GZ(JJ)+SINT*GR(JJ))*AMYU0/A;  % ! AMYU0/A
                    AA(I,sum(NCCN)*2+KNN+K-1+JJ) = AA(I,sum(NCCN)*2+KNN+K-1+JJ)+EE;  %! ### 1
                end
            end
        end
%CC
        fprintf('%s\n','���艻��̉Q�d���� CCS�ɍ�郵');
        for III = 1:CCS
           for I = 1:NCCN(III)
               A = RCCN(III,I);
               B = ZCCN(III,I);
               for K = 1:2:KSN-2
                   [GW,GR,GZ] =  EXTINDC(A,B,RES(K),ZES(K),RES(K+1),ZES(K+1),RES(K+2),...
                   ZES(K+2),NONC,fid99,fid100); % OK
                   for JJ = 1:3
                       EE = GW(JJ)*AMYU0;   %  ! flux*AMYU0
                       AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+KNN+K-1+JJ) = AA(I+sum(NCCN(1:III-1))+NAPB+NFLX,sum(NCCN)*2+KNN+K-1+JJ)+EE; % ! ### 3
                   end
               end 
           end  
       end
       JJJJ = sum(NCCN)+sum(NCCN)+sum(KNN)+sum(KSN);
%
% ??????????????????????????????????????????????????????????????????    
%
         
    end
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

fclose(fid99);
fclose(fid100);

% TOTAL IP(USHIKI)
%C      WRITE(IPR,16) NFLX+NAPB+NCCS
% if (ipconst == 1)
%     fprintf(WAHAHA,'No. of information data = %d\r\n',NFLX+NAPB+sum(NCCN)+1);% +1�͑��d��
%     fprintf(WAHAHA,'Vector FF(I) at the initial stage of iteration\r\n');
%     for I = 1:NFLX+NAPB+sum(NCCN)+1
%         fprintf(WAHAHA,'%d\r\n',FF(I));
%     end
%     fprintf(WAHAHA,'No. of unknowns =   %d   %d   %d   %d\r\n\r\n',JJJJ,sum(NCCN),sum(KNN),sum(KSN));
%     fprintf(WAHAHA,'Matrix AA(I,J)\r\n');
%     for I = 1:NFLX+NAPB+sum(NCCN)+1
%  	    fprintf(WAHAHA,'I=%d\r\n',I);
%         for J = JJJJ
%             fprintf(WAHAHA,'%d\r\n',AA(I,J));
%         end
%     end
% else    
% %C      WRITE(IPR,16) NFLX+NAPB+NCCS
%     fprintf(WAHAHA,'No. of information data = %d\r\n',NFLX+NAPB+sum(NCCN));
%     fprintf(WAHAHA,'Vector FF(I) at the initial stage of iteration\r\n');
%     for I = 1:NFLX+NAPB+sum(NCCN)
%         fprintf(WAHAHA,'%d\r\n',FF(I));
%     end
%     fprintf(WAHAHA,'No. of unknowns =   %d   %d   %d   %d\r\n\r\n',JJJJ,sum(NCCN),KNN,KSN);
%     fprintf(WAHAHA,'Matrix AA(I,J)\r\n');
%     for I = 1:NFLX+NAPB+sum(NCCN)
%  	    fprintf(WAHAHA,'I=%d\r\n',I);
%         for J = JJJJ
%             fprintf(WAHAHA,'%d\r\n',AA(I,J));
%         end
%     end
% end
% ??????????????????????????????????????????????????????????????????    
%
end
FF(1+NAPB:NFLX+NAPB) = FF(1+NAPB:NFLX+NAPB)-GETA+OLDGT;    %! �v�Ċm�F
FF(1+NAPB+NFLX:sum(NCCN)+NAPB+NFLX) = FF(1+NAPB+NFLX:sum(NCCN)+NAPB+NFLX)-GETA+OLDGT;%    ! �v�Ċm�F
end

%% FORM kokomade
%%
%%
%%
%%
%%




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
    fprintf(fid99,'%s\n','/'); 
    fprintf(fid100,'%s\n','/'); 
    fprintf(fid99,'%s\n', '****************************************************');
    fprintf(fid99,'%s\n', '***    In the Subr. INTER **************************');
    fprintf(fid99,'%s\n', '****************************************************');
    fprintf(fid99,'%s\n','/'); 
    fprintf(fid100,'%s\n', '****************************************************');
    fprintf(fid100,'%s\n', '***    In the Subr. INTER **************************');
    fprintf(fid100,'%s\n', '****************************************************');
    fprintf(fid100,'%s\n','/'); 
%
    for I = 1:CCS
        RCCS(I,NCCS(I)+1) = RCCS(I,1);
        ZCCS(I,NCCS(I)+1) = ZCCS(I,1);
    end
%
    DFI(1:sum(NCCN)) = FFOUT(1:sum(NCCN));
    FI(1:sum(NCCN)) = FFOUT(1 + sum(NCCN):sum(NCCN) + sum(NCCN));
    fprintf('%s %d\n', 'GETA in INTER =',GETA);
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
%	         fprintf(WAHAHA,'%s\r\n','Reproducibility of field sensor signals');
% 
	         fid104 = fopen('output/Comparison_TotalFieldSignal.txt','w');%104
	         fid106 = fopen('output/Discrepant_Field_Points.txt','w');%106
             II = 0;
             for I = 1:NAPB
 	             XBBB = XBBR(I)*cos(TET(I+NFLX))+XBBZ(I)*sin(TET(I+NFLX));
%                 fprintf(WAHAHA, '%d Measured = %d Calculated = %d\r\n', I,BSNSR(I),XBBB);
                 fprintf(fid104,'%d %d\n',BSNSR(I),XBBB);
                 EE = 100.0D0*abs((XBBB-BSNSR(I))/BSNSR(I));
                 if (EE > 30.0D0)
                     II = II+1;
                     fprintf(fid106,'%d %d\n', RS(I+NFLX),ZS(I+NFLX));
 %                    fprintf(WAHAHA,'%d %d %d %d %e\r\n', RS(I+NFLX),ZS(I+NFLX),BSNSR(I),XBBB,EE);
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
             fid107 = fopen('output/Discrepant_Flux_Points.txt','w');%107
%           	 fprintf(WAHAHA,'%s\r\n','Reproducibility of flux loop signals');
	         fid105 = fopen('output/Comparison_TotalFluxSignal.txt','w');%105
	         II = 0;
             for I = 1:NFLX
                 fprintf(WAHAHA, '%d Measured = %d Calculated = %d\n', I,FLXLP(I),XPSI(I));
                 fprintf(fid105,'%d %d\n', FLXLP(I),XPSI(I));
                 EE = 100.0D0*abs((XPSI(I)-FLXLP(I))/FLXLP(I));
                 if (EE > 30.0D0)
                     II = II+1;
                     fprintf(fid107,'%d %d\n',RS(I),ZS(I));
                     fprintf(WAHAHA,'%d %d %d %d %d\n', RS(I),ZS(I),FLXLP(I),XPSI(I),EE);
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
   %% �R�C���d������鎥��̉��Z
   [PPSI(1:KCMX),PHIR,PHIZ,PPSIA(1:KCMX),PPSIB(1:KCMX),PIRA,PIRB,PIZA,PIZB,GSTAR,HSTAR,DAG,DBG,DAH,DBH]...
   = STARB(1,AS,BS,RC(1:KCMX),ZC(1:KCMX),RNOR,ZNOR); % OK
   PSIS  = GETA + sum(PPSI(1:KCMX).*ECI(1:KCMX).*RMYU0);
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
    for III = 1:CCS
       for K = 1:NE(III)
           [HW,GW,GR,GZ,HR,HZ] = INTEGS(AS,BS,RCCS(III,2*K-1),ZCCS(III,2*K-1),RCCS(III,2*K),...
           ZCCS(III,2*K),RCCS(III,2*K+1),ZCCS(III,2*K+1)); % OK
           for JJ = 1:3
               KJ2 = 3*(K-1)+JJ + 3*sum(NE(1:III-1));
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
