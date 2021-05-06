function CCS_UTST_matlab(inputfile)
    %% *****************************************************************
    %% *  CAUCHY CONDITION SURFACE METHOD  (quadratic element version) *
    %% *****************************************************************
    
    format long e;
    
    %MXCCS = 20;          % MAX! / NUMBER OF ELEMENTS ON THE CCS
    %MXINT = 10000;       % MAX! / NUMBER OF INTERNAL POINTS
    
    PARAM = loadinputfile(inputfile);
    
    REF =   loadreference(PARAM);

    % [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS] = loadsensordata(PARAM);
    
    [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS] = loadRealsensordata(PARAM);

    % CCS�ʂ������Őݒ�
    for i=1:PARAM.CCS
        PARAM.Z0(i)=    CCS(i);
    end
    
    WALL =  loadwalldata(PARAM);
    
    ExtCOIL = loadcoildata(PARAM);
    
    FFDAT = makeFFdata(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP);

    
    GHR = zeros(1,300); 
    GHZ = zeros(1,300); 
    CCSDAT = makeCCSdata(PARAM, GHR, GHZ);
    
    
    FF = FFDAT;
    FF(end+1:end+sum(CCSDAT.NCCN)) = 0.0;
    if PARAM.IPCONST
        FF(end+1) = -66410*4.0*pi*1.0e-7;
    end
    
    %% *************************************************************
    %% *************************************************************
    %%      iteration Start
    
    DELGE = 0.0;
    OLDEL = 1.0;
    GETA = 0.0;     %! ���ʂ̏����l���[���ɂ���
    OLDGT = GETA;
    OMGA = 1.9;
    
    %fid90 = fopen([PARAM.temporary_file_directory '/GETA.txt'],'w'); %90
    
    if (PARAM.ITSKP > 0)
        ITMX=1;
        fprintf('���ʃT�[�`�� (1) SVD_MT�����݂̂ł��܂�\n');
    else
        fprintf('���ʃT�[�`�� (0) INTER�̔����ł����܂�\n');
    end
    
    % ITEND = 0;
    % EPS1 = 0;
    
    % for IT = 1:ITMX
    %     GETA = GETA + DELGE;
    %     GETA = OMGA*GETA + (1.0D0-OMGA)*OLDGT;
    
    %     fprintf('****************************************************\n');
    %     fprintf('%s%d%s\n','** (IT) Iteration count =',IT,' **');
    %     fprintf('%s%d%d\n','** (GETA/DELGE) =',GETA,DELGE);
    
        %******************************************************************************************************
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %[FC,BR,BZ,PSIFLX,PSIC,AMYU0,AA,FF,GETA,NCCS,fid99,fid100] = FORM(WAHAHA,IT,OLDGT,...
        %GETA,FF,AA,n100,n50,ECI,KCMX,RS,ZS,RC,ZC,ITYPE,TET,NTPB,NNPB,NAPB,NFLX,...
        %RCCN,ZCCN,NCCN,NCCS,RCCS,ZCCS,REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,Nedp,NONC,MXCCS,...
        %RMYU0,NE,CCS,ipconst);%NE uhiski OK
        
        if (PARAM.IPCONST == 1)
            NMAX = SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN)+1;% ushikiip
        else
            NMAX = SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM + sum(CCSDAT.NCCN);% ushikiip
        end
        JMAX = sum(CCSDAT.NCCN) + sum(CCSDAT.NCCN) + sum(WALL.KNN) + sum(WALL.KSN);
        
    %    AA=zeros(length(FF),150);
        AA=zeros(NMAX, JMAX);
    %    AA=zeros(length(FF),sum(CCSDAT.NCCN)*2+WALL.KNN);
        
    %      [FC,BR,BZ,PSIFLX,PSIC,AA,FF] =...
    %      FORM(PARAM,AA,FF,IT,OLDGT,GETA,ExtCOIL,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL);
         [FC,BR,BZ,PSIFLX,PSIC,AA,FF] =...
         FORM(PARAM,AA,FF,ExtCOIL,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL);
     
        %% INOMOTO start
        fluxfactor=10;
    
        FF(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) =    FF(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) * fluxfactor;
        AA(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM,:) =  AA(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM,:) * fluxfactor;
        FC(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) =    FC(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) * fluxfactor;
    
        FLXLP = FFDAT(SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+1:SENSOR_NPRB.NUM+SENSOR_TPRB.NUM+SENSOR_FLXLP.NUM) * fluxfactor;
    
        %% INOMOTO end
    
        %% Solve Matrix Equation
        %% Singular Value Decompoition
    
        fprintf('NCCN/KNN/KSN = %d %d %d\n',sum(CCSDAT.NCCN),sum(WALL.KNN),sum(WALL.KSN));
        
        %  Modified Truncated SVD for Smoothness
        %[C,W,U,V,FFOUT,XBFR,XMT,XGETA,GET] = SVD_MT(ITSKP,IT,AA,FF,NMAX,JMAX,n100,n50,...
        %    0,0.0D0,NAPB,NFLX,NCCN,KNN,KSN,FC,FLXLP,BSNSR,GETA_YN,AUTO,AUTOINPUT);
    %     [C,W,U,V,FFOUT,XBFR,XMT,XGETA,GET] = ...
    %     SVD_MT_matlab(PARAM,PARAM.ITSKP,IT,GETA,AA,FF,FC,0,0.0D0,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL,FLXLP);
        [C,W,U,V,FFOUT,XBFR,XMT] = ...
        SVD_MT_matlab(PARAM,AA,FF,FC,0,0.0D0,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL,FLXLP);
        
        %���m����
        half_norm = sqrt((sum((FFOUT).^2)));
        fprintf('%s%d\r\n','norm of the solution vector = ',half_norm);
    
        %% ***************************
        %% ***************************
    %     if (PARAM.ITSKP > 0)
    %         GETA = XGETA;
    %         break
    %     else
    %     end
        
        %% C***************************
        %% C***************************
    %     [PSI,DELGE,RCCS,ZCCS,XPSI] = INTER(0,GETA,RCCS,ZCCS,CR,CZ,FFOUT,n100,n50,RC,ZC,...
    %     ECIGRP,ECI,KCMX,NAPB,NFLX,FLXLP,BSNSR,RS,ZS,TET,RCCN,ZCCN,...
    %     NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,WAHAHA,...
    %     MXCCS,MXINT,NCCS,NINT,RMYU0,NE,CCS);
        %% C***************************
        %% C***************************
    
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
        
    % end
    
    %  ******************************************************************************************************
    %      iteration End
    %  ******************************************************************************************************
    
    % fprintf('********* Congraturation!!  GETA iteration Converged.\n');
    % fprintf('%s%d\n','****** No. of iterations required = ',ITEND);
    % fprintf('%s%d %d\n','** (GETA/DELGE) = ',GETA,DELGE);
    % fprintf('%s%d\n','** EPS=DABS(DELGE/GETA)= ',EPS1);
    
    %% *************************************************************
    %% *************************************************************
    
    EDDYP(FFOUT,PARAM,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL);
    
    % plot eddy current
    DISF = dlmread([PARAM.output_file_directory '/EddyCurrentProfile.txt']);
    figure('Name','Eddy Current Plofile','NumberTitle','off')
    plot(DISF(:,1),DISF(:,2),'-ko','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',2)
    xlabel({'Distance (m)'});
    ylabel({'Eddy Current Density (MA/m)'});
    axis([0 DISF(end,1) -0.2 0.2])
    set(gca, 'FontSize',14);
    
    axis([0 DISF(end,1) -0.4 0.4])
    set(gca, 'FontSize',14);
    hold on
    JEDDY = dlmread(strcat([PARAM.input_file_directory '/jeddy.txt']));    
    plot(JEDDY(:,1),JEDDY(:,2),'-b')
    
    
    if 0
        % plot sensor position
        VV = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);
        SEN0 = dlmread([PARAM.temporary_file_directory '/SENPOS0.txt']);
        SEN1 = dlmread([PARAM.temporary_file_directory '\SENPOS1.txt']);
        
        figure('Name','Sensor Position','NumberTitle','off')
        subplot(1,2,1);
        plot(VV(:,1),VV(:,2),'-k');
        hold on
        plot(SEN0(:,1),SEN0(:,2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4)
        plot(SEN1(:,1),SEN1(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4)
        title('Sensor Position')
        xlabel('r [m]');
        ylabel('z [m]');
        axis equal
        axis([0 1 -1.2 1.2])
        set(gca, 'FontSize',14);
    
        % plot CV segment
        %COIL = dlmread('output\@UTST_CoilGeom.txt');
        VVMESH = dlmread([PARAM.temporary_file_directory '\VacuumVesselMeshPoints.txt']);
        VVNODE = dlmread([PARAM.temporary_file_directory '\VacuumVesselNodePoints.txt']);
        VVSEG = dlmread([PARAM.temporary_file_directory '\VacuumVesselSegments.txt']);
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
    
    % [PSI,DELGE,RCCS,ZCCS,XPSI] = INTER(1,GETA,RCCS,ZCCS,CR,CZ,FFOUT,n100,n50,RC,ZC,...
    % ECIGRP,ECI,KCMX,NAPB,NFLX,FLXLP,BSNSR,RS,ZS,TET,RCCN,ZCCN,...
    % NCCN,KNE,KNN,REV,ZEV,KSE,KSN,RES,ZES,Nedp,AMYU0,NONC,fid99,fid100,WAHAHA,...
    % MXCCS,MXINT,NCCS,NINT,RMYU0,NE,CCS);
    
    MINR=10;
    MAXR=90;
    MINZ=-100;
    MAXZ=100;
    ICRE=1;
    JCRE=2;
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
    [PSI,DELGE,RCCS,ZCCS,XPSI] = INTER(PARAM,1,GETA,CR,CZ,FFOUT,ExtCOIL,SENSOR_TPRB,SENSOR_NPRB,SENSOR_FLXLP,CCSDAT,WALL,NINT);
    
            CCR = unique(CR);
            CCZ = unique(CZ);
            CCR(1) = [];
            psi=reshape(PSI(1:numel(CCR)*numel(CCZ)),numel(CCZ),numel(CCR));
            figure
            contour(CCR,CCZ,psi,'-k','LevelStep',0.0003); 
            hold on
            contour(REF.R,REF.Z,REF.Flux,'--m','LevelStep',0.0003);  % ���� 
            xlabel({'r (m)'});
            ylabel({'z (m)'});
            axis equal
    
    %% Data positions for 2D plot
    PLOT.R = 0.1:0.01:0.9;
    PLOT.Z = -1:0.02:1;
    
    fclose('all');
    end
    %% Main kokomade