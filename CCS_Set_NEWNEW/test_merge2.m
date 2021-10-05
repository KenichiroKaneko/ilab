var1 = load("vars_beforeLongFor");
var2 = load("vars_afterLongFor");

jt_center = var1.jt_center;
rr = var1.rr;
r = var1.r;
zz = var1.zz;
z = var1.z;
env = var1.env;
psi_virtualj = rr * 0;

lenC = length(jt_center);

% r_center = [0 r(2:1 + lenC) 0];
% jt_center = reshape(jt_center, [1, 1, lenC]);
% jt_center = ones(size(psi_virtualj)) .* jt_center;

% kk = 4 .* rr .* r_center ./ ((rr + r_center + 1e-6).^2 + (zz - z(env.Nz - 1)).^2);
% [K, E] = ellipke(kk);
% psi_virtualj = psi_virtualj + 2 * pi * 4 * pi * 1e-7 .* jt_center ./ (pi * sqrt(kk)) .* sqrt(rr .* r_center) .* ((1 - kk / 2) .* K - E);

% kk = 4 .* rr .* r_center ./ ((rr + r_center + 1e-6).^2 + (zz + z(env.Nz - 1)).^2);
% [K, E] = ellipke(kk);
% psi_virtualj = psi_virtualj + 2 * pi * 4 * pi * 1e-7 .* jt_center ./ (pi * sqrt(kk)) .* sqrt(rr .* r_center) .* ((1 - kk / 2) .* K - E);

% psi_virtualj = sum(psi_virtualj, 3);
% % psi_virtualj(:, 1) = 0;
% % psi_virtualj(:, end) = 0;
% sum(sum(psi_virtualj))

% error('error description')
% const = 2 * pi * 4 * pi * 1e-7;
% tic

% for k = 1:length(jt_center)
%     kk1 = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz - z(env.Nz - 1)).^2);
%     [K1, E1] = ellipke(kk1);
%     kk2 = 4 * rr * r(k + 1) ./ ((rr + r(k + 1) + 1e-6).^2 + (zz + z(env.Nz - 1)).^2);
%     [K2, E2] = ellipke(kk2);

%     % psi_virtualj = psi_virtualj + const * jt_center(k) ./ (pi * sqrt(kk1)) .* sqrt(rr * r(k + 1)) .* ((1 - kk1 / 2) .* K1 - E1);
%     % psi_virtualj = psi_virtualj + const * jt_center(k) ./ (pi * sqrt(kk2)) .* sqrt(rr * r(k + 1)) .* ((1 - kk2 / 2) .* K2 - E2);
%     denominator = ((1 ./ (pi * sqrt(kk1)) .* ((1 - kk1 / 2) .* K1 - E1)) + (1 ./ (pi * sqrt(kk2)) .* ((1 - kk2 / 2) .* K2 - E2)));
%     psi_virtualj = psi_virtualj + const * jt_center(k) .* sqrt(rr * r(k + 1)) .* denominator;
% end

% toc
% sum(sum(psi_virtualj))

psi_virtualj = rr * 0;
const = 2 * pi * 4 * pi * 1e-7;
r_center = r(2:1 + lenC);
tic

for i = 1:25
    % 24こずつにちぎる
    left = 25 * (i - 1) + 1;
    right = left + 23;
    len = 24;

    r_center = reshape(r_center(left:right), [1, 1, len]);
    r_center = ones(size(psi_virtualj)) .* r_center;
    jt_center = reshape(jt_center(left:right), [1, 1, len]);
    jt_center = ones(size(psi_virtualj)) .* jt_center;

    kk1 = 4 * rr .* r_center ./ ((rr + r_center + 1e-6).^2 + (zz - z(env.Nz - 1)).^2);
    [K1, E1] = ellipke(kk1);
    kk2 = 4 * rr .* r_center ./ ((rr + r_center + 1e-6).^2 + (zz + z(env.Nz - 1)).^2);
    [K2, E2] = ellipke(kk2);

    denominator = ((1 ./ (pi * sqrt(kk1)) .* ((1 - kk1 / 2) .* K1 - E1)) + (1 ./ (pi * sqrt(kk2)) .* ((1 - kk2 / 2) .* K2 - E2)));
    psi_virtualj = psi_virtualj + const .* jt_center .* sqrt(rr .* r_center) .* denominator;

end

toc
psi_virtualj = sum(psi_virtualj, 3);
sum(sum(psi_virtualj))
