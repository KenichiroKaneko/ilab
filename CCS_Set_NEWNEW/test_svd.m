% SVDCMPのテスト
vars = load("vars_afterSVDCMP");
format longG
% W = vars.W;
% U = vars.U;
% V = vars.V;

% uu = vars.uu;
% ss = vars.ss;
% vv = vars.vv;

% Ws = size(W)
% Us = size(U)
% Vs = size(V)

% uus = size(uu)
% sss = size(ss)
% vvs = size(vv)
% A = vars.A;
% P = vars.FFOUT;
% F = vars.FC;

% 変数定義１ AP = F という逆問題を解く
% A = [1 2 4 0; 5 -1 1 -1; -1 2 -1 1];
% P = [1; 3; 4; 2];
A = [10 4; 4 8; 2 2];
P = [1.5; 1.25];

[M N] = size(A);

% 変数定義２
M = 10;
N = 2;
A = randi([1 10], M, N)
P = randi([1 10], N, 1);

% 既知の値 Ap = q の q
F = A * P;

% 特異値分解
[U S V] = svd(A)

% 特異値行列の逆行列を計算 S(MxN)→S*(NxM)
S_diag = diag(S);
S_inv = zeros(N, M);

for i = 1:length(S_diag);
    S_inv(i, i) = 1 / S_diag(i);
end

P;
P_solved = V * S_inv * inv(U) * F;
[P P_solved]

error('error description solved')
% 分解→もとのベクトルに戻す
% U * S * V'
% A

% 特異値の小さいやつを消していくぞ
S_diag = diag(S)
S_length = length(S_diag);
KUP = S_length;

for i = 1:S_length

    % 特異値行列の逆行列を計算 S(MxN)→S*(NxM)
    S_inv = zeros(N, M);

    for i = 1:KUP;
        S_inv(i, i) = 1 / S_diag(i);
    end

    % 求めるベクトル
    P_solved = V * S_inv * inv(U) * F;

    % 残差のノルムの計算
    residual(i) = sqrt(sum(sum((A * P_solved - F).^2)));

    % 小さい特異値を減らす
    % S_diag(KUP) = 0;
    KUP = KUP -1;
end

figure()
scatter(1:S_length, residual)
