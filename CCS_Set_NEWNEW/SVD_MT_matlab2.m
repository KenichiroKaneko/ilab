function [C, W, U, V, X, XBFR, XMT] = SVD_MT_matlab2(PARAM, CONFIG, A, B, FC, ...
        LTikh, GAM0, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL, FLXLP, FF)

    [M, N] = size(A);
    NAPB = SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    NFLX = SENSOR_FLXLP.NUM;
    NCCN = CCSDAT.NCCN;
    KNN = WALL.KNN;
    KSN = WALL.KSN;

    % 2021/05/24
    BSNSR = SENSOR_TPRB.TPRB;
    % BSNSR = SENSOR_TPRB.TNPRB;
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
    IPRN = fopen([PARAM.temporary_file_directory '/@SVDCHK.txt'], 'w');

    fprintf(IPRN, '%s\r\n', '   M = NAPB+NFLX+NCCN');
    fprintf(IPRN, '%s%d %d %d %d\r\n', '  = ', M, NAPB, NFLX, sum(NCCN));
    fprintf(IPRN, '%s\r\n', '   N = NCCN+NCCN+KNN+KSN');
    fprintf(IPRN, '%s%d %d %d %d\r\n', '  = ', N, sum(NCCN), sum(NCCN), KNN, KSN);

    fprintf('M/N=%d %d\r\n', M, N);
    fprintf(IPRN, 'M/N=%d %d\r\n', M, N);

    U = A;
    % U = [A zeros(M, M - N)];これ違う

    % 特異値分解 UWV [U S V] = svd(A) に対応するのは W:S特異値1xM, V:V, U:U
    [W, V, U] = SVDCMP(U); % OK
    [uu, ss, vv] = svd(A);
    % [U, S, V] = svd(A);
    % W = diag(S);
    save("vars_afterSVDCMP");
    % error('error description')

    LRSVCHK = 0;

    fprintf(IPRN, '%s\r\n', 'Left SV Vector の check をしません');
    fprintf('%s\r\n', 'Left SV Vector の check をしません');

    LRSVCHK = 1;

    fprintf(IPRN, '%s\r\n', 'Right SV Vector の check をします');
    fprintf('%s\r\n', 'Right SV Vector の check をします');

    %   !  Find maximum singular value
    for K = 1:N
        fprintf(IPRN, '%d\r\n', W(K));
    end

    % ###############################################################################
    % ↓特異値をソート、最大値で割って規格化したものがSVS
    [SVS] = SVSORT_matlab(PARAM, W); % OK

    if CONFIG.ShowFig
        figure('Name', 'Singular Value', 'NumberTitle', 'off')
        semilogy(1:numel(SVS), SVS(1:numel(SVS)), 'ko', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 2)
        % hold on;
        ylabel({'Normalized Singular Value'});
        set(gca, 'FontSize', 14);
    end

    % Xが求める未知数のベクトル1xN (Ap = q の p)
    [C, X, GET] = KUPCHK_matlab(PARAM, A, B, U, V, W, NAPB, NFLX, sum(NCCN)); % OK ushiki

    % save('vars_KUPCHK');
    % error('error description kupchk');

    ITRNC = PARAM.ITRNC;
    IDCN = PARAM.IDCN; % 1

    % KUP0 = PARAM.KUP0;
    % 2021/05/10
    % KUP0 = curvature(SVS);
    % scatter(KUP0, log(SVS(KUP0)), 'ro')
    % hold off;

    % L-curve法 2021/05/17
    % KUP0 = LCURVE(PARAM, CONFIG, A, ss, vv, uu, X, FC);
    % KUP0 = 30;
    % KUP0 = 50;
    % KUP0 = 62;
    % KUP0 = N;
    KUP0 = 65;

    % 2021/06/11 各要素の相対誤差を計算１
    CalcMRE(PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, A, ss, vv, uu, X, FC, FF);

    fprintf(IPRN, '%s %d %s\r\n', 'You truncate SVs smaller than', KUP0, '-th SV');
    fprintf('%s %d %s\r\n', 'You truncate SVs smaller than', KUP0, '-th SV');

    fid209 = fopen([PARAM.temporary_file_directory '/@Determined_Gap.txt'], 'w'); % 209
    CKUP = KUP0 + KUP0 + 1;
    CKUP = CKUP / 2.0D0;
    fprintf(fid209, '%d %d\r\n', CKUP, W(N));
    fprintf('%s\r\n', 'OK?');

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
    fprintf(IPRN, 'Condition Number = %d\r\n', CONDNO);
    fprintf('Condition Number = %d\r\n', CONDNO);

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
    %        !  Zero the "small" singular values
    % 打切りしている↓
    W(KUP0 + 1:N) = 0.0D0;
    KUP = KUP0;

    fprintf('%s %d\r\n', '***  KUP=', KUP);
    fprintf('%s %d\r\n', '***  KUP0=', KUP0);

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
    for K = 1:N
        fprintf(IPRN, '%d\r\n', W(K));
    end

    fprintf('Improved Condition Number = %d\r\n', CONDNO);
    CONDNO = SMAX00 / SMIN;
    fprintf('Newly Defined Improved Cond.No. = %d\r\n', CONDNO);
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
    % 2021/06/11 各要素の相対誤差を計算２
    CalcMRE(PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, A, ss, vv, uu, X, FC);
    % save('vars_afterSVBKSB')
    % error('error description svbksb')

    fprintf(IPRN, '%s\r\n', '    Solution vector is:');

    for K = 1:N
        fprintf(IPRN, '%d\r\n', X(K));
    end

    fprintf(IPRN, '%s\r\n', '     Original right-hand side vector:');

    for I = 1:M
        fprintf(IPRN, '%d\r\n', C(K));
    end

    fprintf(IPRN, '%s\r\n', '     Result of (matrix)*(sol''n vector):');
    %
    C = zeros(1, M);
    C(1:M) = A(1:M, 1:N) * X(1:N)';
    %
    %***************************
    %***************************
    %        !↓ここ、下駄の件で手直しがいる。
    RESD(1:M) = (C(1:M) - (B(1:M) - XGETA)) .* and(1:M > NAPB, 1:M <= (NAPB + NFLX)) ...
        + (C(1:M) - B(1:M)) .* or(1:M <= NAPB, 1:M > (NAPB + NFLX));
    %***************************
    %***************************
    for K = 1:M
        fprintf(IPRN, '%d\r\n', C(K));
    end

    fprintf(IPRN, '%s\r\n', ' ');
    fid101 = fopen([PARAM.temporary_file_directory '/Comparison_FieldSignal.txt'], 'w'); %101
    fid102 = fopen([PARAM.temporary_file_directory '/Comparison_FluxSignal.txt'], 'w'); %102
    fid201 = fopen([PARAM.temporary_file_directory '/Comparison_TotField.txt'], 'w'); %201
    fid202 = fopen([PARAM.temporary_file_directory '/Comparison_TotFlux.txt'], 'w'); %202
    fid103 = fopen([PARAM.temporary_file_directory '/Line45deg.txt'], 'w'); %103
    VAL = 1.0D03;
    fprintf(fid103, '%d %d\r\n', -VAL, -VAL);
    fprintf(fid103, '%d %d\r\n', +VAL, +VAL);
    fprintf(IPRN, '%s\r\n', 'The given & the reconstructed signals');
    %
    for K = 1:M

        if (K > NAPB && K <= (NAPB + NFLX))
            fprintf(IPRN, '%d %d %d %d\r\n', K, B(K) - XGETA, C(K), RESD(K));
        else
            fprintf(IPRN, '%d %d %d %d\r\n', K, B(K), C(K), RESD(K));
        end

        if (K <= NAPB)
            fprintf(fid101, '%d %d\r\n', B(K), C(K));
            fprintf(fid201, '%d %d\r\n', BSNSR(K), C(K) + FC(K));
        else

            if (K <= NAPB + NFLX)
                fprintf(fid102, '%d %d\r\n', B(K) - XGETA, C(K)); % ! B(K),C(K)共に下駄が差し引かれている。
                fprintf(fid202, '%d %d\r\n', FLXLP(K - NAPB), C(K) + FC(K) + GETA); % ! 測定値 FLXLP() には下駄が含まれる。
            end

        end

        if or(K == NAPB, K == (NFLX + NAPB))
            fprintf(IPRN, '%s\r\n', '/');
        end

    end

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

    fid84 = fopen([PARAM.temporary_file_directory '/SmallSVs.txt'], 'w');

    if (KUP < N)

        for I = KUP + 1:N
            fprintf(fid84, '%d %d\r\n', I, SVS(I));
        end

    else

        for I = 1:N
            fprintf(fid84, '%d %d\r\n', I, 1.0D0 - 20);
        end

    end

    fclose(IPRN);
    fclose(fid84);
    fclose(fid101);
    fclose(fid102);
    fclose(fid103);
    fclose(fid201);
    fclose(fid202);
    fclose(fid209);
end
