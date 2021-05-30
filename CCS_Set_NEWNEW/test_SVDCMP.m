format longG
A = randi([1 100], 4, 3)
% A = [A; 0 0 0 0]

[W, V, U] = SVDCMP(A)

[u s v] = svd(A)

W = sort(W, 'descend')

% SVDCMPのテスト
vars = load("vars_afterSVDCMP");
format longG
A = vars.A;
W = vars.W;
U = vars.U;
V = vars.V;

uu = vars.uu;
ss = vars.ss;
vv = vars.vv;

Ws = size(W)
Us = size(U)
Vs = size(V)

uus = size(uu)
sss = size(ss)
vvs = size(vv)

W = diag(W);

sqrt(sum(sum((U * W * V - A).^2)))
W(90, 87) = 0;
sqrt(sum(sum((uu * ss * vv - A).^2)))
