function makeRealSenposCCSdata(psi, Bz, Br)

    f0 = readmatrix("sensorCoordinate0.txt");
    f1 = readmatrix("sensorCoordinate1.txt");

    Nz = 2033;
    Nr = 800;
    zmin = -9.985000000000001e-01;
    zmax = 9.985000000000001e-01;
    rmin = 1.081500000000000e-01;
    rmax = 0.8880;
    delr = 9.747920133111480e-04;
    delz = 9.827755905511811e-04;

    zz = linspace(zmin, zmax, Nz);
    rr = linspace(rmin, rmax, Nr);

    psiOutside = load("psiOutside.mat").A;

    padding = zeros(198, 2033);
    psi602 = [psi; padding];
    size(psi602)

    psi = psi602 + psiOutside';

    % 磁束
    r = f0(:, 1);
    z = f0(:, 2);
    psi_CCS = psi(r, z);
    fp = fopen( "Sensor_B.txt", "w");
    fprintf(fp, 'r[m]\tz[m]\tpsi[Wb]\tBz[T]\tBr[T]\n');
    for i = 1:length(r)
        fprintf(fp, '%f\t%f\t%f\t0\t0\n', rr(r(i)), zz(z(i)), psi_CCS(i));
    end
    fclose(fp)
    % sensor_B = [r, z, psi_CCS];
    % writematrix(sensor_B, "test_sensor_B.txt");

    % 磁場
    r = f1(:, 1);
    z = f1(:, 2);
    Bz_CCS = Bz(r, z);
    Br_CCS = Br(r, z);
    fp = fopen( "Sensor_Flux.txt", "w");
    fprintf(fp, 'r[m]\tz[m]\tpsi[Wb]\tBz[T]\tBr[T]\n');
    for i = 1:length(r)
        fprintf(fp, '%f\t%f\t0\t%f\t%f\n', rr(i), zz(i), Bz_CCS(i), Br_CCS(i));
    end
    fclose(fp)
    % sensor_flux = [r, z, Bz_CCS, Br_CCS];
    % writematrix(sensor_flux, "test_sensor_flux.txt");

end
