%% FORM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% FORM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% FORM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% FORM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% AA*X=FF
%function [FC,BR,BZ,PSIFLX,PSIC,AMYU0,AA,FF,GETA,NCCS,fid99,fid100]= FORM(WAHAHA,...
%IT,OLDGT,GETA,FF,AA,n100,n50,ECI,KCMX,RS,ZS,RC,ZC,ITYPE,TET,NTPB,NNPB,...
%NAPB,NFLX,RCCN,ZCCN,NCCN,NCCS,RCCS,ZCCS,REV,ZEV,KNE,KNN,RES,ZES,KSE,KSN,Nedp,...
%NONC,MXCCS,RMYU0,NE,CCS,ipconst)

function [FC,BR,BZ,PSIFLX,PSIC,AA,FF] =...
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
