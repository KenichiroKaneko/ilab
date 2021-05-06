%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS] = loadsensordata(PARAM)

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
        
        if ((R< 0.57+0.0015 && R>0.27-0.0015)||(Z<0.239 && R>0.67)||(R<0.689 && R>0.629-0.0015))
        else
            chnum = chnum+1;
            
            SENSOR_TPRB.R(chnum*2-1) =       R;
            SENSOR_TPRB.R(chnum*2) =         R;
            SENSOR_TPRB.Z(chnum*2-1) =       Z;
            SENSOR_TPRB.Z(chnum*2) =         -Z;
            SENSOR_TPRB.TET(chnum*2-1) =     atan2(BZ,BR);
            SENSOR_TPRB.TET(chnum*2) =       atan2(BZ,-BR);
            SENSOR_TPRB.TPRB(chnum*2-1) =    sqrt(BR^2 + BZ^2);
            SENSOR_TPRB.TPRB(chnum*2) =      sqrt(BR^2 + BZ^2);
            SENSOR_TPRB.ITYPE(chnum*2-1) =   1;
            SENSOR_TPRB.ITYPE(chnum * 2) = 1;
            RR(chnum) = R;
            ZZ(chnum) = Z;
            B(chnum) = sqrt(BR^2 + BZ^2);
        end
    end

    SENSOR_TPRB.NUM = chnum*2;
    disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);

    figure()
    plot(SENSOR_TPRB.R, SENSOR_TPRB.Z, 'o');
    title("Sensor B pos");
    % error('fin')
    % インボード側の磁束のプロットと、CCS面の決定2021年4月18日
    figure()
    RR = [flip(RR) RR ];
    ZZ = [-flip(ZZ) ZZ];
    B = [flip(B) B ];
    RR0index = RR < 0.15;
    B = B(RR0index)
    ZZ = ZZ(RR0index);
    RR
    ZZ
    plot(B, ZZ);
    lmin = islocalmin(B);
    lmax = islocalmax(B);
    hold on
    plot(B(lmin), ZZ(lmin), 'r*');
    hold on
    plot(B(lmax), ZZ(lmax), 'b*');
    title("r=0")
    lmin = ZZ(lmin);
    matrix = [B(lmax)' ZZ(lmax)'];
    matrix = sortrows(matrix, 'descend');
    % error("test")

    if length(ZZ(lmax)) > 1

        for i = 1:2
            CCS(i) = matrix(i, 2);
        end

    else
        size(matrix)
        CCS(1) = matrix(1, 2);
    end

    % ここまで

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


    R = 0;
    Z = 0;

    for i = 1:sensornum_B
        R(i) = str2double(sensordata_B{6 + (i - 1) * 5});
        Z(i) = str2double(sensordata_B{7 + (i - 1) * 5});
    end

    figure()
    plot(R, Z, 'o');
    title("sensor flux pos")
    R = 0;
    Z = 0;

    chnum = 0;

    for i=1:sensornum_Flux
        R = str2double(sensordata_Flux{6+(i-1)*5});
        Z = str2double(sensordata_Flux{7+(i-1)*5});
        PSI = str2double(sensordata_Flux{8+(i-1)*5});
        BZ = str2double(sensordata_Flux{9+(i-1)*5});
        BR = str2double(sensordata_Flux{10+(i-1)*5});

        chnum = chnum+1;

        SENSOR_FLXLP.R(chnum*2-1) =       R;
        SENSOR_FLXLP.R(chnum*2) =         R;
        SENSOR_FLXLP.Z(chnum*2-1) =       Z;
        SENSOR_FLXLP.Z(chnum*2) =         -Z;
        SENSOR_FLXLP.FLXLP(chnum*2-1) =   PSI;
        SENSOR_FLXLP.FLXLP(chnum*2) =     PSI;
        SENSOR_FLXLP.TET(chnum*2-1) =     0.0D0;
	    SENSOR_FLXLP.TET(chnum*2) =       0.0D0;
	    SENSOR_FLXLP.ITYPE(chnum*2-1) =   0;
	    SENSOR_FLXLP.ITYPE(chnum*2) =     0;
    end
    
    SENSOR_FLXLP.NUM=chnum*2;
    disp(['Number of FLXLP = ' num2str(SENSOR_FLXLP.NUM)]);

    figure()
    plot(SENSOR_FLXLP.R, SENSOR_FLXLP.Z, 'o');
    title("sensor flux loop")
    

%     %% Write files for CCS
%     fp = fopen([PARAM.temporary_file_directory '/SENPOS0.txt'],'w'); % 110
%     for i=1:SENSOR_FLXLP.NUM
%         fprintf(fp,'%d %d\n',SENSOR_FLXLP.R(i),SENSOR_FLXLP.Z(i));
%     end
%     fclose(fp);
%     
%     fp = fopen([PARAM.temporary_file_directory '/SENPOS1.txt'],'w'); % 111
%     for i=1:SENSOR_TPRB.NUM
%         fprintf(fp,'%d %d\n',SENSOR_TPRB.R(i),SENSOR_TPRB.Z(i));
%     end
%     fclose(fp);
%     
%     fp = fopen([PARAM.temporary_file_directory '/SENPOS2.txt'],'w'); % 112
%     for i=1:SENSOR_NPRB.NUM
%         fprintf(fp,'%d %d\n',SENSOR_NPRB.R(i),SENSOR_NPRB.Z(i));
%     end
%     fclose(fp);
% 
    if PARAM.IUTST==5
        SENSOR_FLXLP.FLXLP=SENSOR_FLXLP.FLXLP-0.0042459;
    end
%     
%     %%  *************************************************************************
%     %%     Generation of CCS input data
%     %%  *************************************************************************
% 
%     fprintf('Generation of CCS input data ** START ***\n');
%     fp = fopen([PARAM.temporary_file_directory '/CCSinput_UTST(temp).txt'],'w');
%     
%     fprintf(fp,'%s\n','*');
%     fprintf(fp,'%s\n','*** CCS Test input for UTST generated in PreUTST ***'); 
%     fprintf(fp,'%s\n','** NTPB/NNPB/NFLX=(No. of T-/N-probes/Flux-Loops) **');
%     fprintf(fp,'   %d     %d     %d\n',SENSOR_TPRB.NUM,SENSOR_NPRB.NUM,SENSOR_FLXLP.NUM); 
%     fprintf(fp,'%s\n','**  GETA (SSURF)');
%     fprintf(fp,'  %s\n','0.0000E+00');
%     fprintf(fp,'%s\n','* T-Probe');
%     for i = 1:SENSOR_TPRB.NUM
%         fprintf(fp,' %d\n',SENSOR_TPRB.TPRB(i));
%     end
%     fprintf(fp,'%s\n','* N-Probe');
%     fprintf(fp,'%s\n','* Flux-Loop');
%     for i = 1:SENSOR_FLXLP.NUM
%         fprintf(fp,' %d\n',SENSOR_FLXLP.FLXLP(i));
%     end
%     fprintf(fp,'%s\n','****** MINR * MAXR * MINZ * MAXZ ****');
%     fprintf(fp,'%s\n','10   90  -100  100');
%     fprintf(fp,'%s\n','*********');
%     fprintf(fp,'%s\n','* ---�R�C���d���f�[�^�̕���---�P��[kA]');
%     fprintf(fp,'%s\n','* EF');
%     fprintf(fp,'%s\n','* PF#1');
%     fprintf(fp,'%s\n','* PF#2');
%     fprintf(fp,'%s\n','* PF#3');
%     fprintf(fp,'%s\n','* PF#4');
% 
%     fprintf(fp,'   %s\n','-0.28');
%     fprintf(fp,'   %s\n','0.0');
%     fprintf(fp,'   %s\n','0.0');
%     fprintf(fp,'   %s\n','0.0');
%     fprintf(fp,'   %s\n','0.0');
%     
%     fclose(fp);
%     fprintf('Generation of CCS input data ** END ***\n');
end
%% loadsensordata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
