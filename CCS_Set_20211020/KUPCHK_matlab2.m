%function [C,X,GET] = KUPCHK_matlab(PARAM,A,B,U,V,W,M,N,MP,NP,NAPB,NFLX,NCCN)
function [C, X, GET] = KUPCHK_matlab(PARAM, A, B, U, V, W, NAPB, NFLX, NCCN)
    [M, N] = size(A);
    C = zeros(1, N); %MP
    WTMP = zeros(1, 200);
    AIC0 = zeros(1, 200);
    AIC1 = zeros(1, 200);
    GET = zeros(1, N); %NP
    %
    %  ****  AX=B  ****
    %  Input==>  A(M,N): System Matrix,
    %            B(M): Inhomogeneous Term Vector
    %  Output==> X(N): Unknown Vector to be Solved
    %   M=Number of data points or equations,   N=Number of unknowns
    %   MP=Max. capacity of M,                  NP=Max. capacity of N
    % *****************************************************************
    %
    %
    %
    AI0MIN = 1.0D10;
    AI1MIN = 1.0D10;

    for KUP = 1:N
        WTMP(1:N) = W(1:N);

        if (KUP ~= N)
            WTMP(KUP + 1:N) = 0.0D0;
        end

        %********** GETA iteration (START) *****************
        ITMX = 1 + PARAM.GETA_YN * 999;
        %ITMX = 1000;
        DELGE = 0.0D0;
        OLDEL = 1.0D0;
        GETA = 0.0D0;
        OLDGT = GETA;
        OMGA = 1.9D0;
        %           OMGA = 1.D0;
        for IT = 1:ITMX
            GETA = GETA + DELGE;
            GETA = OMGA * GETA + (1.0D0 - OMGA) * OLDGT;
            C(1:M) = B(1:M);
            C(NAPB + 1:NAPB + NFLX) = C(NAPB + 1:NAPB + NFLX) - GETA;
            %
            [X] = SVBKSB(U, WTMP, V, C); %OK
            K = 1:M;
            J = 1:N;
            C(K) = A(K, J) * X(J)';
            E(K) = (C(K) - B(K) + GETA) .* and(K > NAPB, K <= (NFLX + NAPB)) + (C(K) - B(K)) .* or(K <= NAPB, K > (NFLX + NAPB)); % ! GETAを引くのはFlux loopのみ !!
            DELGE = sum((B(K) - C(K) - GETA) .* and(K > NAPB, K <= (NFLX + NAPB))); %  ! DELGE評価はFlux loopのみ
            BNSNALL = sum(E(K) .* E(K));
            BNSN_B = sum(E(K) .* E(K) .* (K <= NAPB));
            BNSN_P = sum(E(K) .* E(K) .* and(K > NAPB, K <= (NFLX + NAPB)));
            BNSN_C = sum(E(K) .* E(K) .* (K > (NFLX + NAPB)));
            CNT0 = NFLX;
            DELGE = DELGE / CNT0;
            EPSDEL = abs((DELGE - OLDEL) / OLDEL);

            OLDEL = DELGE;
            OLDGT = GETA;
            EPS = abs(DELGE / GETA);
            ITEND = IT;

            if (EPS < 1.0E-5)
                break
            end

        end

        GET(KUP) = GETA;
    end

end
