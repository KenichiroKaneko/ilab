function [KCMX,RC,ZC,ECI] = COIL_UTST_2(INP,WAHAHA)
%    
% ========================================================================
%     coil data(position, makisuu).
%    
% ------------------------------------------------------------------------
%     RCOIL  : R position.
%     ZCOIL  : Z position.
%     CCOIL  : coil vector.
%     MCOIL  : coil makisu.
%    
% =========================================================================
% NAC = 5;
% EFの実形状を考慮
NAC = 10;
ECIGRP = zeros(NAC);
RCOIL = zeros(2,NAC); ZCOIL = zeros(2,NAC); CCOIL = zeros(2,NAC);
%
% MCOIL  : coil makisu (Turns).
% EF位置をZ方向に分散させる
%   coil counter -->   EF1  EF2  EF3  EF4  EF5  EF6    PF#1    PF#2    PF#3    PF#4
MCOIL              = [ 200, 200, 200, 200, 200, 200,   8,      3,      8,      3];  
%
%     RCOIL ( R-POSTION OF COIL ) --------------------------------------
%               bottom    upper
%...EF1-COIL
RCOIL(:,1) =  [ 0.80D0,   0.80D0];
%...EF2-COIL
RCOIL(:,2) =  [ 0.80D0,   0.80D0];
%...EF3-COIL
RCOIL(:,3) =  [ 0.80D0,   0.80D0];
%...EF4-COIL
RCOIL(:,4) =  [ 0.80D0,   0.80D0];
%...EF5-COIL
RCOIL(:,5) =  [ 0.80D0,   0.80D0];
%...EF6-COIL
RCOIL(:,6) =  [ 0.80D0,   0.80D0];
%
%...PF#1-COIL
RCOIL(:,7) =  [ 0.20D0,   0.20D0];
%
%...PF#2-COIL
RCOIL(:,8) =  [ 0.685D0,  0.685D0];
%
%...PF#3-COIL
RCOIL(:,9) =  [ 0.75D0,   0.75D0];
%
%...PF#4-COIL
RCOIL(:,10) =  [ 0.685D0,  0.685D0];
%
%     ZCOIL ( Z-POSITION OF COIL ) ------------------------------------
%               bottom    upper
%以前はEFコイルのZ位置を1.07mに集約して計算していた→EF実構造を考慮し、6分割して計算
%...EF1-COIL
ZCOIL(:,1) =  [ -1.13D0,   1.13D0];
%...EF2-COIL
ZCOIL(:,2) =  [ -1.11D0,   1.11D0];
%...EF3-COIL
ZCOIL(:,3) =  [ -1.09D0,   1.09D0];
%...EF4-COIL
ZCOIL(:,4) =  [ -1.05D0,   1.05D0];
%...EF5-COIL
ZCOIL(:,5) =  [ -1.03D0,   1.03D0];
%...EF6-COIL
ZCOIL(:,6) =  [ -1.010D0,   1.01D0];
%
%...PF#1-COIL
ZCOIL(:,7) =  [-1.10D0,   1.10D0];
%
%...PF#2-COIL
ZCOIL(:,8) =  [-0.80D0,   0.80D0];
%
%...PF#3-COIL
ZCOIL(:,9) =  [-0.675D0,  0.675D0];
%
%...PF#4-COIL
ZCOIL(:,10) =  [-0.50D0,   0.50D0];
%
%     CCOIL ( COIL VECTER ) ----------------------------------
%
%...EF-COIL
CCOIL(:,1) =  [ 1.0D0,    1.0D0];
%...EF-COIL
CCOIL(:,2) =  [ 1.0D0,    1.0D0];
%...EF-COIL
CCOIL(:,3) =  [ 1.0D0,    1.0D0];
%...EF-COIL
CCOIL(:,4) =  [ 1.0D0,    1.0D0];
%...EF-COIL
CCOIL(:,5) =  [ 1.0D0,    1.0D0];
%...EF-COIL
CCOIL(:,6) =  [ 1.0D0,    1.0D0];
%
%...PF#1-COIL
CCOIL(:,7) =  [ 1.0D0,    1.0D0];
%
%...PF#2-COIL
CCOIL(:,8) =  [ 1.0D0,    1.0D0];
%
%...PF#3-COIL
CCOIL(:,9) =  [ 1.0D0,    1.0D0];
%
%...PF#4-COIL
CCOIL(:,10) =  [ 1.0D0,    1.0D0];
%
% *******************************************************************

for I = 1:NAC
    HEAD = textscan(INP,'%1s',1); 
    while (HEAD{1}{1} == '*')
        textscan(INP,'%s',1,'Delimiter','\n'); %806
        HEAD = textscan(INP,'%1s',1);
        
    end
    temp = textscan(INP,'%s',1,'Delimiter','\n');
    %コイル電流値[kA]
    disp(HEAD{1}{1})
    disp(temp{1}{1})
%     
    
    ECIGRP(I) = str2double(strcat(HEAD{1}{1},temp{1}{1}));
%     ECIGRP(I) = str2double(strcat(HEAD{1}{1},temp{1}{1}));
        
    
%%cd	ECIGRP(I)=ECIGRP(I)/(2.0D0*PAI)
%     ECIGRP(I) = -ECIGRP(I)*1000.0D0; % ECIGRP=Current of each coil group [kA]-->[A]  これ JT-60U
    ECIGRP(I) = ECIGRP(I)*1000.0D0; % -がつくのはおかしい？
    disp(ECIGRP)
    
%%CC      ECIGRP(I)=+ECIGRP(I)*1000.0D0 ! ECIGRP=Current of each coil group [kA]-->[A]
end %                       !*** 座標系は（R,-φ,Z）? ***prof. inomoto の作成した磁束プロファイルの座標系に合わせるため？→ 前の正解磁束プロファイルはプラズマ電流負にしていた
%++++++++++++++++++++++++++++++++++++++++++++++++++++
fprintf(WAHAHA,'%s\r\n','Coil data has been successfully installed!');   
%++++++++++++++++++++++++++++++++++++++++++++++++++++
I = 0;
RC  = zeros(sum(MCOIL));
ZC  = zeros(sum(MCOIL));
ECI = zeros(sum(MCOIL));
for N = 1:NAC % 100
    if (abs(ECIGRP(N)) < 1.0D-10)
        continue
    end
	MM = MCOIL(N);
    for KM = 1:2
      for M = 1:MM % 110
%c      DO 110 M=1,1
%Caution!!   All turns are assumed to be at the same position.
%電流源の位置を設定（現段階では電流源が同位置に存在）
	      I = I+1;
	      RC(I)  = RCOIL(KM,N);
	      ZC(I)  = ZCOIL(KM,N);
	      ECI(I) = CCOIL(KM,N)*ECIGRP(N);
%Caution!!      
          IPRNT = 1;
          if (IPRNT <= 0)
              continue
          end
	      fprintf(WAHAHA,'%d Type=%d %d RS/ZS= %d %d  ECI= %d\r\n',I,N,M,RC(I),ZC(I),ECI(I));
          if (M == MM)
             fprintf(WAHAHA,' \r\n');
          end
      end % 110
   end
end %100
KCMX = I;
    