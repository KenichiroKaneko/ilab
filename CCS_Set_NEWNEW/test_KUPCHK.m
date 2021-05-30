vars = load("vars_KUPCHK");
vars2 = load("vars_afterForm");
vars3 = load("vars_afterSVBKSB");

format longG
M = 90;
N = 87;

X = vars.X;
A = vars.A; % 状態方程式
U = vars.U; % USVのU
V = vars.V; % USVのV
W = vars.W; % USVのS、対角行列の対角成分、特異値
FF = vars2.FF; % センサーの値

[u S v] = svd(A);

S_diag = diag(S);
S_inv = zeros(N, M);

for i = 1:length(S_diag);
    S_inv(i, i) = 1 / S_diag(i);
end

P;
P_solved = v * S_inv * inv(u) * FF';

% [P' P_solved X']
sum(sum((A * X' - FF).^2))
sum(sum((A * P_solved - FF).^2))

X2 = vars3.X;

sum(sum((A * X2' - FF).^2))
