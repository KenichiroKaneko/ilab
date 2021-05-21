function [CXK,CYK,SOLK,MSK] = FNGRPH_UTST(IUTST,eddy_time,SensorVer)
format long e
NCX=101;
NCY=201;
if (IUTST == 0)
elseif (IUTST == 1)
SOLK = dlmread('@UTST_FluxProfile.txt');
elseif (IUTST == 2)
SOLK = dlmread('@UTST_FluxProfile_#2.txt');
elseif (IUTST == 3)
SOLK = dlmread('@UTST_FluxProfile_#3.txt');
elseif (IUTST == 4)
SOLK = dlmread('@UTST_FluxProfile_#4.txt'); %! 0731
elseif (IUTST == 5)
%vacuum�̂��߂ɊԎ؂�
% SOLK = dlmread('@UTST_FluxProfile_#12.txt'); % 
% SOLK = dlmread('@UTST_FluxProfile_#5.txt'); %! 0731
SOLK = dlmread('@UTST_FluxProfile_#5_1.txt'); %! kondo GETA処理を要しない
NCX=602;
NCY=2033;
elseif (IUTST == 6)
SOLK = dlmread('@UTST_FluxProfile_#6.txt'); % 1201
% SOLK = dlmread('@UTST_FluxProfile_#6_1.txt'); %! kondo GETA処理を要しない（見直し）
NCX=602;
NCY=2033;
elseif (IUTST == 7)
SOLK = dlmread('@UTST_FluxProfile_#7.txt'); % 1201
NCX=602;
NCY=2033;
elseif (IUTST == 8)
SOLK = dlmread(strcat('@UTST_FluxProfile_#8_',num2str(eddy_time),'.txt')); % 1201
NCX=602;
NCY=2033;
elseif (IUTST == 9)
SOLK = dlmread('@UTST_FluxProfile_#9.txt'); % 
NCX=602;
NCY=2033;
% elseif (IUTST == 10)
% SOLK = dlmread('@UTST_FluxProfile_#10_m.txt'); % 
% NCX=602;
% NCY=2033;
elseif (IUTST == 10)
SOLK = dlmread('@UTST_FluxProfile_#10_m2.txt'); % 
NCX=602;
NCY=2033;
elseif (IUTST == 11)
SOLK = dlmread('@UTST_FluxProfile_#10.txt'); % 
NCX=602;
NCY=2033;
elseif (IUTST == 12)
SOLK = dlmread('@UTST_FluxProfile_#12.txt'); % 
NCX=604;
NCY=2033;
end
CXK  = zeros(1,NCX);
CYK  = zeros(1,NCY);
MSK  = zeros(NCX,NCY,'int8');
%
IRMX=NCX;
IZMX=NCY;
%
%*****
for IZ=1:IZMX
    for IR=1:IRMX
        SOLK(IZ,IR) = SOLK(IZ,IR)/(2.0*pi); % ���C�ʊ֐��̊֌W
    end
end
%
ZMIN=-0.9985 ;
ZMAX=0.9985;
RMIN=0.108150  ;
RMAX=0.6940;
CRMX=IRMX-1;
CZMX=IZMX-1;
DELZ=(ZMAX-ZMIN)/CZMX;
DELR=(RMAX-RMIN)/CRMX;
for IZ = 1:IZMX
    CIZM1 = IZ-1;
    Z = ZMIN+DELZ*CIZM1;
    for IR=1:IRMX
        CIRM1 = IR-1;
        R = RMIN+DELR*CIRM1; %R�̗��U���@0.108150�{(0.6940-0.108150)/(602-1)*CIRM1
	    CXK(IR) = R;
	    CYK(IZ) = Z;
	    MSK(IR,IZ) = 0;
        if and( abs(Z) > 0.285, R > 0.5985) 
            SOLK(IZ,IR) = 0.0; % ���̕ǂ̊O�͕`���Ȃ�
        end
    end
end
end      