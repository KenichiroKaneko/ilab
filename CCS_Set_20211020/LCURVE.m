% Lcurve法で打切り特異値分解の打切り項数を決定する

function KUP = LCURVE(PARAM, CONFIG, A, W, V, U, X, FC)
    % A：特異値分解される行列
    % W：特異値の行列
    % V：[U S V] = svd(A)のV
    % U：[U S V] = svd(A)のU
    % X：求める未知変数のベクトル
    % FC：Ap = q の q
    % M：既知の数、N：未知の数

    [M N] = size(A);

    FC = FC';
    W = diag(W); % 対角成分のみの行列
    % W = sort(W, "descend");
    W_length = length(W);
    truncateNum = W_length;

    for i = 1:W_length

        % 特異値行列の逆行列を計算 S(MxN)→S*(NxM)
        W_inv = zeros(N, M);

        for i = 1:truncateNum;
            W_inv(i, i) = 1 / W(i);
        end

        % 求めるベクトル
        P_solved = V * W_inv * inv(U) * FC;

        % 残差のノルムの計算
        residual(i) = norm(A * P_solved - FC).^2;

        % 求めたベクトルのノルム
        P_norm(i) = norm(P_solved).^2;

        % 小さい特異値を減らす
        truncateNum = truncateNum -1;
    end

    if 1 %CONFIG.ShowFig

        figure()
        loglog(residual, P_norm, 's')
        % ax = gca;
        % ax.XDir = 'reverse';
        title("L-curve");
        xlabel("残差 ||Ap^* - q||")
        ylabel("||p^*||")

        % i = 10;
        % text(residual(i), (P_norm(i)), ['\leftarrow', num2str(i)])

        for i = 5:5:W_length
            text(residual(i), (P_norm(i)), ['\leftarrow', num2str(i)])
        end

    end

    % error('error description')
    % save("vars_inLcurve");

    KUP = CURVATURE(residual, P_norm)

    % error('error description lcurve')

end
