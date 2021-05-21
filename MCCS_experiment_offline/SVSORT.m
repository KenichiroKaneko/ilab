function [SVS] = SVSORT(N,SVO,NMX)
SVS = zeros(NMX);
fid13 = fopen('@SV_SORT.txt','w');%13
fid14 = fopen('@SV_STEEP.txt','w');%14
SVMX = -1.0D10;	
for K = 1:N
	SVS(K) = SVO(K);
    if (SVS(K) > SVMX)
	    SVMX = SVS(K);
    else
    end
end
fprintf('%s %d\r\n','SVMX =',SVMX);
for J = 2:N
	A = SVS(J);
    for I = J-1:-1:1
        if (SVS(I) >= A)
            break
        end
	    SVS(I+1) = SVS(I);
        I = 0;
    end
    SVS(I+1) = A;
end

SVS(1:N) = SVS(1:N)/SVMX;
%SVS(1:N) = SVS(1:N);
for K = 1:N
    fprintf(fid13,'%d %d\r\n',K,SVS(K));
end

for K = 1:N-1
	EE = SVS(K+1)/SVS(K);
	STEEP = log(EE);
    fprintf(fid14,'%d %d\r\n',K,STEEP);
end