function makeRealSenposCCSdata(psi, Bz, Br)

    f0 = readmatrix("sensorCoordinate0.txt");
    f1 = readmatrix("sensorCoordinate1.txt");

    psiOutside = load("psiOutside.mat").A;

    padding = zeros(198, 2033);
    psi602 = [psi; padding];
    size(psi602)

    psi = psi602 + psiOutside';

    % 磁束
    r = f0(:, 1);
    z = f0(:, 2);
    psi_CCS = psi(r, z);
    sensor_B = [r, z, psi_CCS];
    writematrix(sensor_B, "test_sensor_B.txt");

    % 磁場
    r = f1(:, 1);
    z = f1(:, 2);
    Bz_CCS = Bz(r, z);
    Br_CCS = Br(r, z);
    sensor_flux = [r, z, Bz_CCS, Br_CCS];
    writematrix(sensor_flux, "test_sensor_flux.txt");

end
