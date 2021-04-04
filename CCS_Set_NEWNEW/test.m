 [M,N]=size(A)
 UU=zeros(M,M);
 UU(1:M,1:N)=U;
 WW=diag(W);
% WW2=zeros(N,M);
% WW2(1:N,1:N)=WW;

 WW3=zeros(M,N);
 WW3(1:N,1:N)=WW;

 Wt=W(1:KUP);
 WWt=diag(Wt);
 WW3t=zeros(M,N);
 WW3t(1:KUP,1:KUP)=WWt;

% AA2=V*WW2*UU;
 AA3=UU*WW3*V';
 
% norm(A-AA2')
 norm(A-AA3)