%  [M,N]=size(A)
%  UU=zeros(M,M);
%  UU(1:M,1:N)=U;
%  WW=diag(W);
% % WW2=zeros(N,M);
% % WW2(1:N,1:N)=WW;

%  WW3=zeros(M,N);
%  WW3(1:N,1:N)=WW;

%  Wt=W(1:KUP);
%  WWt=diag(Wt);
%  WW3t=zeros(M,N);
%  WW3t(1:KUP,1:KUP)=WWt;

% % AA2=V*WW2*UU;
%  AA3=UU*WW3*V';

% % norm(A-AA2')
%  norm(A-AA3)

% filenames = ["z2033_r602" "z1000_r602" "z0960_r602" "z0920_r602" ...
    %              "z0880_r602" "z0840_r602" "z0800_r602" "z0760_r602" ...
    %              "z0720_r602" "z0680_r602" "z0640_r602" "z0600_r602" ...
    %              "z0560_r602" "z0520_r602" "z0480_r602" "z0440_r602" ...
    %              "z0400_r602" "z0360_r602" "z0320_r602" "z0280_r602" ...
    %              "z0240_r602" "z0200_r602" "z0160_r602" "z0120_r602" ...
    %              "z0080_r602"];
% folder_name = "newfolder2021";

% for i = 1:length(filenames)

%     str = extractAfter(filenames(i), "z");
%     str = extractBefore(str, "_");
%     folder_name = "UTST_numel_" + str;

%     if not(exist(folder_name, 'dir'))
%         mkdir(folder_name)
%     end

% end

GSplot_CCS_for_finemesh_merge2("z2033_r602")
