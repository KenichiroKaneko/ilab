% Adopted_num�iAUTOINPUT(31)�j���Œ�


function [C,W,U,V,X,XBFR,XMT,GETA,GET] = SVD_MT_190623(ITSKP,IT,A,B,M,N,MP,NP,LTikh,...
GAM0,NAPB,NFLX,NCCN,KNN,KSN,FC,FLXLP,BSNSR,GETA_YN,AUTO,AUTOINPUT, Adopted_num)
% B=FF
U = zeros(MP,NP);
C = zeros(1,MP);	
RESD = zeros(1,2100);
XBFR = zeros(1,NP);
XMT = zeros(1,NP);
KUP0 = 0; %ushiki
%
%  ****  AX=B  ****
%  Input==>  A(M,N): System Matrix, 
%            B(M): Inhomogeneous Term Vector
%  Output==> X(N): Unknown Vector to be Solved
%   M=Number of data points or equations,   N=Number of unknowns
%   MP=Max. capacity of M,                  NP=Max. capacity of N
% *****************************************************************
%
IPRN =   fopen('@SVDCHK.txt','w');
fprintf('%s%d\r\n','IT=',IT);
if (IT > 1) %1001
else
    fprintf(IPRN,'%s\r\n','   M = NAPB+NFLX+NCCN');
    fprintf(IPRN,'%s%d %d %d %d\r\n','  = ',M,NAPB,NFLX,sum(NCCN));
    fprintf(IPRN,'%s\r\n','   N = NCCN+NCCN+KNN+KSN');
    fprintf(IPRN,'%s%d %d %d %d\r\n','  = ',N,sum(NCCN),sum(NCCN),KNN,KSN);
    fprintf('MP/NP/M/N=%d %d %d %d\r\n' ,MP,NP,M,N);
    fprintf(IPRN,'MP/NP/M/N=%d %d %d %d\r\n' ,MP,NP,M,N);
    if (N > NP)
	    fprintf('%s\r\n', 'Lack of Dimension (N>NP) in Subr. SVD_MT');
	    return%    CALL EXIT(0)
    else
    end
    if (M > MP)
	    fprintf('%s\r\n', 'Lack of Dimension (M>MP) in Subr. SVD_MT');
	    return %CALL EXIT(0)
    else
    end
%
%    !   Copy A into U
    U(1:M,1:N) = A(1:M,1:N);
%    !  Decompose matrix A
%   W:Singuler value
    [W,V,U] = SVDCMP(U,M,N,MP,NP); % OK
% 
    LRSVCHK = 0;
    if (LRSVCHK > 0)     
        fprintf(IPRN,'%s\r\n','Left SV Vector �� check �����܂�');
        fprintf('%s\r\n','Left SV Vector �� check �����܂�');
%        !Check Left Singular vectors  (Columns of U(I,J))
%        for J = 1:M 
%            WRITE(ClNo,"(I4)") J+1000
%            fname1="VECTOR_U"//ClNo(2:4)//".txt"
%            OPEN(191,file=fname1)
%            for I = 1:M 
%                WRITE(191,*) I,U(I,J)
%            end
%c      PAUSE 'OK?'
%        end
    else
        fprintf(IPRN,'%s\r\n','Left SV Vector �� check �����܂���');
        fprintf('%s\r\n','Left SV Vector �� check �����܂���');
    end      
%
    LRSVCHK = 1;
    if (LRSVCHK > 0)     
        fprintf(IPRN,'%s\r\n','Right SV Vector �� check �����܂�');
        fprintf('%s\r\n','Right SV Vector �� check �����܂�');
%        !Check Right Singular vector (Columns of V(I,J))
%         for J = 1:N 
%             WRITE(ClNo,"(I4)") J+1000
%             fname2="VECTOR_V"//ClNo(2:4)//".txt"
%             OPEN(192,file=fname2)
%             for I = 1:N 
%                 WRITE(192,*) I,V(I,J)
%             end
% %c      PAUSE 'OK?'
%         end   
    else
        fprintf(IPRN,'%s\r\n','Right SV Vector �� check �����܂���');
        fprintf('%s\r\n','Right SV Vector �� check �����܂���');
    end
%   !  Find maximum singular value
    for K = 1:N
        fprintf(IPRN,'%d\r\n',W(K));
    end
% ###############################################################################
    [SVS] = SVSORT(N,W,NP); % OK
%     figure('Name','Singular Value','NumberTitle','off')
%     semilogy(1:numel(SVS), SVS(1:numel(SVS)), 'ko','MarkerEdgeColor','k',...
%                                       'MarkerFaceColor','k','MarkerSize',2)
%     ylabel({'Normalized Singular Value'});
%     set(gca, 'FontSize',14);
%    fprintf('%s\r\n','   ');
%    fprintf('%s\r\n','***** Now check Singular values & AIC *****');
    if (IT == 1)                                                 
        [C,X,GET] = KUPCHK(A,B,U,V,W,M,N,MP,NP,NAPB,NFLX,sum(NCCN));  % OK ushiki    
    end
%CCC
    if (AUTO == 0)
        prompt = 'Truncate SMALL singular values (SVs)? (Yes/No)=(1/0)\n';
        ITRNC = input(prompt);
    else
        ITRNC = AUTOINPUT(28);
    end
    if (ITRNC <= 0)      
        COND0 = 1.0D30;%  ! Non-Truncated SVD
        fprintf('%s\r\n', 'SVs Not Truncated');
    else
        if (AUTO == 0)
            prompt = 'Input, Cond.No.(0) or No. of Untruncated SVs (1)?\n';
            IDCN = input(prompt);
        else
            IDCN = AUTOINPUT(29);
        end
        if (IDCN == 0)
            prompt = 'Then, Input your desired Cond. No. (6000.?)\n';
            COND0 = input(prompt);
            fprintf(IPRN,'%s %d\r\n','Small SVs are Truncated so that Cond.No.  < ',COND0);
            fprintf(IPRN,'%s %d\r\n','Small SVs are Truncated so that Cond.No.  < ',COND0);
        else
            if (AUTOINPUT(30) == 1)
                AUTO = 0;
            end
            if or(AUTO == 0, AUTOINPUT(30) == 0)           
                prompt = 'Then, Input No. of Untruncated Singular Values (98?)\n';
                KU
                P0 = input(prompt);
            else
                KUP0 = Adopted_num;
            end
            fprintf(IPRN,'%s %d %s\r\n','You truncate SVs smaller than', KUP0,'-th SV');
            fprintf('%s %d %s\r\n','You truncate SVs smaller than', KUP0,'-th SV');
        end
        fid209 = fopen('@Determined_Gap.txt','w'); % 209
        CKUP = KUP0+KUP0+1;
        CKUP = CKUP/2.0D0;
%            fprintf(fid209,'%d %d\r\n' ,CKUP,W(KUP0));ushiki
        fprintf(fid209,'%d %d\r\n', CKUP,W(N));
        fprintf('%s\r\n','OK?');
        if (AUTO == 0)
%             pause;
        end
    end
% ###############################################################################
%
    SMAX = -1.0E20;
	SMIN = +1.0E20;
    WMAX = 0.0;
    for K = 1:N
        if (W(K) > WMAX)
            WMAX = W(K);
        end
        if (W(K) > SMAX)
            SMAX = W(K);
        end
        if (W(K) < SMIN);
            SMIN = W(K);
        end
    end
    CONDNO = SMAX/SMIN;
    SMAX00 = SMAX;
	fprintf(IPRN,'Condition Number = %d\r\n',CONDNO);
	fprintf('Condition Number = %d\r\n',CONDNO);
%    !  Define "small"
%!      WMIN=WMAX*(1.0E-6)
    if (ITRNC <= 0)
        KUP = N;
    else
        if (LTikh == 0)
            if (IDCN == 0)    
                COND = COND0;
            else
                COND = WMAX/W(KUP0);
            end
        else
            COND = 1.0d30;
        end
%cccccccccccccc
        WMIN = WMAX/COND;
%        !  Zero the "small" singular values
        if (IDCN == 0) 
            KUP = 0;
            for K = 1:N
                if (W(K) < WMIN)
	                W(K)=0.0D0;
                else
	                KUP = KUP+1;
                end
            end
        else
          	W(KUP0+1:N) = 0.0D0;
            KUP = KUP0;
        end 
%        WRITE(12,*) '***  KUP=',KUP
        fprintf('%s %d\r\n','***  KUP=',KUP);
        fprintf('%s %d\r\n','***  KUP0=',KUP0);
        
        if (GETA_YN == 1)
            GETA = GET(KUP);
        else 
            GETA = 0;
        end
%
%cccccccccccccccccccccccccccccccccccccccccccccccc
 %       ! DAMPING!
        GAMMA = 0.0D0.*(LTikh == 0) + GAM0.*(LTikh ~= 0);
%        ! #2
        SMAX = -1.0E20;
	    SMIN = +1.0E20;
        for K = 1:N
            if (W(K) == 0.0D0)          %GOTO 444
                continue
            else
                W(K) = (W(K)^2+GAMMA)/W(K);
                if (W(K) > SMAX)
                    SMAX = W(K);
                elseif (W(K) < SMIN)
                    SMIN=W(K);
                end
            end
        end %444
        CONDNO = SMAX/SMIN;
%cccccccccccccccccccccccccccccccccccccccccccccccc
%c
        for K = 1:N
           fprintf(IPRN,'%d\r\n', W(K));
        end
        fprintf('Improved Condition Number = %d\r\n',CONDNO);
        CONDNO = SMAX00/SMIN;
        fprintf('Newly Defined Improved Cond.No. = %d\r\n',CONDNO);
%
%
% **********************************************************************      
%  Backsubstitute for each right-hand side vector
% **********************************************************************      
%
    end
end
%    !�������A���ʂ̌��Ŏ蒼��������B   
XGETA = GETA.*(ITSKP > 0);
C(1:M) = (B(1:M)-XGETA).*and(1:M > NAPB, 1:M <= (NAPB+NFLX))...
	      + B(1:M).*or(1:M <= NAPB, 1:M > (NAPB+NFLX));
%***************************
%***************************
%
[X] = SVBKSB(U,W,V,M,N,MP,NP,C); % OK
%      
fprintf(IPRN,'%s\r\n','    Solution vector is:');
for K = 1:N
    fprintf(IPRN,'%d\r\n',X(K));
end
fprintf(IPRN,'%s\r\n','     Original right-hand side vector:');
for I = 1:M
    fprintf(IPRN,'%d\r\n',C(K));
end
fprintf(IPRN,'%s\r\n','     Result of (matrix)*(sol''n vector):');
%
C = zeros(1,M);
C(1:M) = A(1:M,1:N)*X(1:N)';
%   
%***************************
%***************************
%        !�������A���ʂ̌��Ŏ蒼��������B   
RESD(1:M) = (C(1:M)-(B(1:M)-XGETA)).*and(1:M > NAPB, 1:M <= (NAPB+NFLX))...
          + (C(1:M)-B(1:M)).*or(1:M <= NAPB, 1:M > (NAPB+NFLX));
%***************************
%***************************
for K = 1:M
    fprintf(IPRN,'%d\r\n',C(K));
end
fprintf(IPRN,'%s\r\n', ' ');
	fid101 = fopen('Comparison_FieldSignal.txt','w');%101
	fid102 = fopen('Comparison_FluxSignal.txt','w');%102
	fid201 = fopen('Comparison_TotField.txt','w');%201
	fid202 = fopen('Comparison_TotFlux.txt','w');%202
	fid103 = fopen('Line45deg.txt','w');%103
	VAL = 1.0D03;
	fprintf(fid103,'%d %d\r\n', -VAL,-VAL);
	fprintf(fid103,'%d %d\r\n', +VAL,+VAL);
	fprintf(IPRN,'%s\r\n','The given & the reconstructed signals');
%
    for K = 1:M
        if (K > NAPB && K <= (NAPB+NFLX))
	        fprintf(IPRN,'%d %d %d %d\r\n',K,B(K)-XGETA,C(K),RESD(K));
        else
	        fprintf(IPRN,'%d %d %d %d\r\n',K,B(K),C(K),RESD(K));
        end
    if (K  <= NAPB)
	    fprintf(fid101,'%d %d\r\n',B(K),C(K)); %B�F�V�O�i���x�N�g���AC�F�č\���l�i�Z���T�ʒu�ɂ�����j
	    fprintf(fid201,'%d %d\r\n',BSNSR(K),C(K)+FC(K));  %BSNSR�FFFDAT�AC�F�č\���l�i�Z���T�ʒu�ɂ�����j�AFC�F�R�C���d����^
    else
        if(K <= NAPB+NFLX)
	        fprintf(fid102,'%d %d\r\n', B(K)-XGETA,C(K)); % ! B(K),C(K)���ɉ��ʂ�����������Ă���B�@B�F�V�O�i���x�N�g���AC�F�č\���l�i�Z���T�ʒu�ɂ�����j
	        fprintf(fid202,'%d %d\r\n',FLXLP(K-NAPB),C(K)+FC(K)+GETA);% ! ����l FLXLP() �ɂ͉��ʂ��܂܂��B�@FLXLP�F�t���b�N�X���[�v�V�O�i���AC+FC�{GETA�F�����ȍč\���l�i�R�C���d����^�A���ʂ��܂ށj
        end
    end
    if or(K == NAPB, K == (NFLX+NAPB))
        fprintf(IPRN,'%s\r\n','/');
    end
    end
    if (LTikh == 0)
        XBFR(1:N) = X(1:N);
% ********************************************	  
% ********************************************	  
%        !  Modified Truncated Singular Value Decomposition by Hansen
        if (IT > 1)
            return
        end
        if (AUTO == 0)
            prompt = 'Ordinary TSVD? / Modified TSVD? =(0/1)\n';
            MTS = input(prompt);
        else
            MTS = AUTOINPUT(32);
        end
        if (MTS > 0)
           if (ITRNC > 0)
%         ! KUP���A����`�ϐ��ɒ��ӂ��Č�������̂��ƁB         
               [XMT,X] = MTSVD(V,X,NP,N,KUP,sum(NCCN),KNN,KSN); % OK
           else
           end
        else
        end
% ********************************************	  
% ********************************************	  
    else
    end
    
      fid84 = fopen('SmallSVs.txt','w'); 
      if (KUP < N)
          for I = KUP+1:N
              fprintf(fid84,'%d %d\r\n', I,SVS(I));
          end
      else
          for I = 1:N
              fprintf(fid84,'%d %d\r\n', I,1.0D-20);
          end
      end
