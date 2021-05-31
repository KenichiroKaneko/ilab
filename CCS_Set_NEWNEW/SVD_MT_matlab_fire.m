function [C, W, U, V, X, XBFR, XMT] = SVD_MT_matlab2(PARAM, A, B, FC, ...
        LTikh, GAM0, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, FLXLP, superKUP)

    [M, N] = size(A);
    NAPB = SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    NFLX = SENSOR_FLXLP.NUM;
    NCCN = CCSDAT.NCCN;
    KNN = WALL.KNN;
    KSN = WALL.KSN;

    % 2021/05/24
    BSNSR = SENSOR_TPRB.TPRB;
    GETA_YN = PARAM.GETA_YN;

    AUTO = 1;

    U = zeros(M, N);
    C = zeros(1, M);
    RESD = zeros(1, 2100);
    XBFR = zeros(1, N);
    XMT = zeros(1, N);
    KUP0 = 0;

    %  ****  AX=B  ****
    %  Input==>  A(M,N): System Matrix,
    %            B(M): Inhomogeneous Term Vector
    %  Output==> X(N): Unknown Vector to be Solved
    %   M=Number of data points or equations,   N=Number of unknowns
    %   MP=Max. capacity of M,                  NP=Max. capacity of N
    % *****************************************************************
    %

    U = A;

    % 特異値分解 UWV [U S V] = svd(A) に対応するのは W:S特異値1xM, V:V, U:U
    [W, V, U] = SVDCMP(U); % OK
    [uu, ss, vv] = svd(A);

    % ###############################################################################
    % ↓特異値をソート、最大値で割って規格化したものがSVS
    [SVS] = SVSORT_matlab(PARAM, W); % OK

    % Xが求める未知数のベクトル1xN (Ap = q の p)
    [C, X, GET] = KUPCHK_matlab2(PARAM, A, B, U, V, W, NAPB, NFLX, sum(NCCN)); % OK ushiki

    ITRNC = PARAM.ITRNC;
    IDCN = PARAM.IDCN; % 1

    % L-curve法 2021/05/17
    % KUP0 = LCURVE(PARAM, A, ss, vv, uu, X, FC);
    KUP0 = superKUP;

    CKUP = KUP0 + KUP0 + 1;
    CKUP = CKUP / 2.0D0;
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

        if (W(K) < SMIN)
            SMIN = W(K);
        end

    end

    CONDNO = SMAX / SMIN;
    SMAX00 = SMAX;

    if (LTikh == 0) % true

        if (IDCN == 0) % false
            COND = COND0;
        else
            COND = WMAX / W(KUP0); % これが実行される
        end

    else
        COND = 1.0d30;
    end

    WMIN = WMAX / COND;
    % 打切りしている↓
    W(KUP0 + 1:N) = 0.0D0;
    KUP = KUP0;

    GETA = 0;

    %
    %cccccccccccccccccccccccccccccccccccccccccccccccc
    %       ! DAMPING!
    GAMMA = 0.0D0 .* (LTikh == 0) + GAM0 .* (LTikh ~= 0);
    SMAX = -1.0E20;
    SMIN = +1.0E20;

    for K = 1:N

        if (W(K) == 0.0D0) %GOTO 444
            continue
        else
            W(K) = (W(K)^2 + GAMMA) / W(K);

            if (W(K) > SMAX)
                SMAX = W(K);
            elseif (W(K) < SMIN)
                SMIN = W(K);
            end

        end

    end %444

    CONDNO = SMAX / SMIN;
    %cccccccccccccccccccccccccccccccccccccccccccccccc
    %c
    CONDNO = SMAX00 / SMIN;
    %
    %
    % **********************************************************************
    %  Backsubstitute for each right-hand side vector
    % **********************************************************************
    %

    %    !↓ここ、下駄の件で手直しがいる。
    %XGETA = GETA.*(ITSKP > 0);
    XGETA = GETA;
    C(1:M) = (B(1:M) - XGETA) .* and(1:M > NAPB, 1:M <= (NAPB + NFLX)) ...
        + B(1:M) .* or(1:M <= NAPB, 1:M > (NAPB + NFLX));
    %***************************
    %***************************
    %
    [X] = SVBKSB(U, W, V, C); % OK

    % save('vars_afterSVBKSB')
    % error('error description svbksb')

    %
    C = zeros(1, M);
    C(1:M) = A(1:M, 1:N) * X(1:N)';
    %
    %
    if (LTikh == 0)
        XBFR(1:N) = X(1:N);
        % ********************************************
        %        !  Modified Truncated Singular Value Decomposition by Hansen
        %         if (IT > 1)
        %             return
        %         end
        if (AUTO == 0)
            prompt = 'Ordinary TSVD? / Modified TSVD? =(0/1)\n';
            MTS = input(prompt);
        else
            MTS = PARAM.MTS;
        end

        if (MTS > 0)

            if (ITRNC > 0)
                %         ! KUP等、未定義変数に注意して後日訂正のこと。
                % MTSVD : Modified Trancated SVD
                [XMT, X] = MTSVD(V, X, NP, N, KUP, sum(NCCN), KNN, KSN);
            end

        end

    end

end
