function [X] = SVBKSB(U,W,V,M,N,MP,NP,B)
NMAX = 4900; %  !12��5�� NNEQ�ƘA��
%CC
X = zeros(1,N);
TMP = zeros(1,N);%NMAX);
if (NMAX >= N)
    TMP(1:N) = (U(1:M,1:N)'*B(1:M)').*(W(1:N)' ~= 0)./W(1:N)'; % 
    TMP(isnan(TMP)) = 0;
    X(1:N) = V(1:N,1:N)*TMP(1:N)';
    
%     V(1:N,1:N) * W(1:N) * U(1:M,1:N)'
    
else   
    fprintf('%s%d %s%d\r\n','NMAX=',NMAX,'   N=',N);
    fprintf('%s\r\n','Lack of Dimension (NMAX<N) in Subr. SVBKSB');
end%
end