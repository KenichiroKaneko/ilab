function[Bz, Psi] = EF_calc(r, z, r_size, z_size, V)

%%%与えられた座標（r, z）におけるEF磁場、磁束を計算

% close all


% r = [0.11, 0.16, 0.25, 0.32, 0.39, 0.46, 0.53, 0.60, 0.67];
% z = [0.23, 0.18, 0.12, 0.06, 0.00, -0.06, -0.12, -0.18, -0.23];

R_EF = 855 * 1e-3;
z_EF = 1.05;
Turn_EF = 200;
mu0 = 4 * pi * 1e-7;
% V = 120;
I = -(0.849 * (1.19 * V - 5.32) - 5.56);

Bz = zeros(z_size, r_size);
Psi = zeros(z_size, r_size);
for i = 1 : z_size    
    alpha_u = sqrt((R_EF + r).^2 + (z(i) - z_EF).^2);
    alpha_l = sqrt((R_EF + r).^2 + (z(i) + z_EF).^2);
    k2_u = (4 * R_EF * r) ./ alpha_u.^2;
    k_u = sqrt(k2_u);
    k2_l = (4 * R_EF * r) ./ alpha_l.^2;
    k_l = sqrt(k2_l);
    m_u = k2_u;
    m_l = k2_l;
    [K_u, E_u] = ellipke(m_u);
    [K_l, E_l] = ellipke(m_l);
    Bz_upper = ((mu0 * Turn_EF * I) ./ (2 * pi)) .* (1 ./ alpha_u) .* (K_u + ((R_EF^2 - r.^2 - (z(i) - z_EF).^2) ./ ((R_EF - r).^2 + (z(i) - z_EF).^2)) .* E_u);
    Bz_lower = ((mu0 * Turn_EF * I) ./ (2 * pi)) .* (1 ./ alpha_l) .* (K_l + ((R_EF^2 - r.^2 - (z(i) + z_EF).^2) ./ ((R_EF - r).^2 + (z(i) + z_EF).^2)) .* E_l);
    Bz(i, :) = Bz_upper + Bz_lower;
    
    %磁束について
    %ベクトルポテンシャルのΦ成分
    A_phai_u = (mu0 * Turn_EF * I)  ./ (pi) .* sqrt(R_EF ./ r) .* (1 ./ k_u) .* ((1 - k_u.^2 ./ 2) .* K_u - E_u);
    A_phai_l = (mu0 * Turn_EF * I)  ./ (pi) .* sqrt(R_EF ./ r) .* (1 ./ k_l) .* ((1 - k_l.^2 ./ 2) .* K_l - E_l);
    
    psi_u = 2 * pi * r .* A_phai_u;
    psi_l = 2 * pi * r .* A_phai_l;
    Psi(i, :) = psi_u + psi_l;
    
%     figure('Name','Flux EF','NumberTitle','off')
%     plot(z, Psi(i, :))
%     pause
end



% Bz = zeros(9, 9);
% for j = 1 : 9
%     for i = 1 : 9
%         r(i);
%         z(j);
%         alpha_u = sqrt((R_EF + r(i))^2 + (z(j) - z_EF)^2);
%         alpha_l = sqrt((R_EF + r(i))^2 + (z(j) + z_EF)^2);
%         k_u = (4 * R_EF * r(i)) / alpha_u^2;
%         k_l = (4 * R_EF * r(i)) / alpha_l^2;
%         m_u = k_u;
%         m_l = k_l;
%         K_u = ellipticCK(m_u);
%         E_u = ellipticCE(m_u);
%         K_l = ellipticCK(m_l);
%         E_l = ellipticCE(m_l);
%         Bz_upper = ((mu0 * Turn_EF * I) / (2 * pi)) * (1 / alpha_u) * (K_u + ((R_EF^2 - r(i)^2 - (z(j) - z_EF)^2) / ((R_EF - r(i))^2 + (z(j) - z_EF)^2)) * E_u);
%         Bz_lower = ((mu0 * Turn_EF * I) / (2 * pi)) * (1 / alpha_l) * (K_l + ((R_EF^2 - r(i)^2 - (z(j) + z_EF)^2) / ((R_EF - r(i))^2 + (z(j) + z_EF)^2)) * E_l);
%         Bz(j, i) = Bz_upper + Bz_lower;
%     end
% %     plot(r, Bz(j, :))
% %     pause
% end
% 
% for i = 1: 9
%     plot(r, Bz(i, :))
%     hold on
% end

% plot(r(1), Bz(1, :))

