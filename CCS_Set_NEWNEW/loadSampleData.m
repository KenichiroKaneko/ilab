%% loadsensordata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R, REF] = loadSampleData(PARAM, CONFIG)

    % CCSセンサー位置を推定する用（R位置は使っていないのでPARAMの値を使用）
    CCS_R = 0;
    CCS_Z = 0;

    % matファイルでの読み込み
    % GSplot_CCS_for_finemesh_merge2()でmerged.matファイルを作成した
    vars = load(PARAM.input_file_directory + "/merged.mat");
    env = vars.env3c;

    % プラズマ無しの容器内の磁束
    psi602 = zeros(env.Nz, env.Nr_orig);
    [env, psi602] = cal_psi(env, psi602);
    psi602 = psi602.';

    % プラズマ無しの容器外も含めた磁束
    env.Nr = 800;
    env.rmax = env.rmin + env.delr * (env.Nr - 1);
    psi800 = zeros(env.Nz, env.Nr);
    [env, psi800] = cal_psi(env, psi800);
    psi800 = psi800.';

    % プラズマ有りの容器内の磁束
    psi_orig = vars.psi;

    % プラズマ有りの容器外も含めた磁束
    padding = zeros(env.Nr - env.Nr_orig, env.Nz);
    psi602 = [psi602; padding];
    psi_orig = [psi_orig; padding];
    psi = psi_orig - psi602 + psi800;
    % psi = psi800;

    % 磁場
    br = vars.Br;
    bz = vars.Bz;
    z = linspace(env.zmin, env.zmax, env.Nz);
    r = linspace(env.rmin, env.rmax, env.Nr);

    REF.Flux = psi';
    REF.R = r;
    REF.Z = z;

    %% 実際のセンサー配置でデータを読み込む
    sensorPosB = fileread("./CCS_temporary/CCS_MP_sensor_position.txt");
    sensorPosB = strsplit(sensorPosB, {'\n', '\t', '\r'});
    sensorPosFL = fileread("./CCS_temporary/CCS_FLXLP_sensor_position_3.txt");
    sensorPosFL = strsplit(sensorPosFL, {'\n', '\t', '\r'});
    SENSOR_TPRB.NUM = (length(sensorPosB) - 5) / 5;
    SENSOR_FLXLP.NUM = (length(sensorPosFL) - 5) / 5;
    % センサー配置をテキストデータから取得
    SENSOR_TPRB.R = str2double(sensorPosB(6:5:5 * SENSOR_TPRB.NUM + 5));
    SENSOR_TPRB.Z = str2double(sensorPosB(7:5:5 * SENSOR_TPRB.NUM + 5));
    SENSOR_FLXLP.R = str2double(sensorPosFL(6:5:5 * SENSOR_FLXLP.NUM + 5));
    SENSOR_FLXLP.Z = str2double(sensorPosFL(7:5:5 * SENSOR_FLXLP.NUM + 5));
    % 磁場や磁束のデータをmatファイルから取得
    % センサー位置の最小値でインデックスを探索している
    for i = 1:SENSOR_TPRB.NUM
        [m, I] = min(abs(r - SENSOR_TPRB.R(i)));
        indexB_R(i) = I;
        [m, I] = min(abs(z -SENSOR_TPRB.Z(i)));
        indexB_Z(i) = I;
    end

    for i = 1:SENSOR_FLXLP.NUM
        [m, I] = min(abs(r - SENSOR_FLXLP.R(i)));
        indexFL_R(i) = I;
        [m, I] = min(abs(z -SENSOR_FLXLP.Z(i)));
        indexFL_Z(i) = I;
    end

    SENSOR_TPRB.R = r(indexB_R);
    SENSOR_TPRB.Z = z(indexB_Z);
    SENSOR_FLXLP.R = r(indexFL_R);
    SENSOR_FLXLP.Z = z(indexFL_Z);

    for i = 1:SENSOR_TPRB.NUM
        BR(i) = br(indexB_R(i), indexB_Z(i));
        BZ(i) = bz(indexB_R(i), indexB_Z(i));
    end

    for i = 1:SENSOR_FLXLP.NUM
        PSI(i) = psi(indexFL_R(i), indexFL_Z(i));
    end

    SENSOR_TPRB.TET = atan2(BZ, BR);
    SENSOR_TPRB.TPRB = sqrt(BR.^2 + BZ.^2);
    SENSOR_TPRB.ITYPE = 1; % 使っていない

    %% No NPRB
    SENSOR_NPRB.NUM = 0;
    SENSOR_NPRB.R = [];
    SENSOR_NPRB.Z = [];
    SENSOR_NPRB.TET = [];
    SENSOR_NPRB.NPRB = [];
    SENSOR_NPRB.ITYPE = [];
    RR = SENSOR_TPRB.R;
    ZZ = SENSOR_TPRB.Z;
    B = sqrt(BR.^2 + BZ.^2);
    SENSOR_TPRB.TNPRB = sqrt(BR.^2 + BZ.^2);
    disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);

    if CONFIG.DevideFlux
        SENSOR_FLXLP.FLXLP = PSI ./ (pi .* SENSOR_FLXLP.R .* SENSOR_FLXLP.R);
    else
        SENSOR_FLXLP.FLXLP = PSI;
    end

    SENSOR_FLXLP.TET = zeros(1, SENSOR_FLXLP.NUM);
    SENSOR_FLXLP.ITYPE = zeros(1, SENSOR_FLXLP.NUM);
    disp(['Number of FLXLP = ' num2str(SENSOR_FLXLP.NUM)]);

    % CCSセンター位置の決定
    if CONFIG.DetermineCCSZPos
        RR0index = RR == min(RR);
        B = B(RR0index); ZZ = ZZ(RR0index);
        ZZ = ZZ'; B = B';
        f = fit(ZZ, B, 'smoothingspline', 'SmoothingParam', 1);
        x = fnzeros(fnder(f.p));
        x = unique(x(:));
        y = fnval(f.p, x);
        matrix = [x y];
        matrix = sortrows(matrix, 2, 'descend');
        lmax = islocalmax(B);

        if CONFIG.ShowFig
            figure()
            hold on
            plot(ZZ, B);
            fnplt(f.p);
            plot(x, y, 'o');
            xlim([-1 1]);
            plot(ZZ(lmax), B(lmax), 'b*');
            view(90, 90);
        end

        if length(ZZ(lmax)) > 1

            for i = 1:2
                CCS_Z(i) = matrix(i, 1);
            end

        else
            CCS_Z(1) = matrix(1, 1);
        end

    else

        for i = 1:PARAM.CCS
            CCS_Z(i) = PARAM.Z0(i);
        end

    end

    % CCS面の自動決定R方向（使っていない）
    if CONFIG.DetermineCCSRPos
        CCS_R = CalcPlasmaCenter(PARAM, CONFIG, CCS_Z);
        % CCS_R = CCS_R * 0.85;
        % CCS_R = 0;

    else

        for i = 1:PARAM.CCS
            CCS_R(i) = PARAM.R0(i);
        end

    end

end

function [param, psi] = cal_psi(param, psi)

    Mu = 4 * pi * 1.0e-7; % 真空の透磁率
    mtrxz = (0:1:param.Nz - 1)' .* param.delz + param.zmin;
    mtrxr = (0:1:param.Nr - 1) .* param.delr + param.rmin;

    for k = 1:param.ncoil
        % (z*rの大きさの配列を作成)
        tmp1 = (param.coil_z(1, k) - mtrxz).^2 * ones(1, param.Nr) + ones(param.Nz, 1) * (param.coil_r(1, k) + mtrxr).^2;
        kk = (4 .* ones(param.Nz, 1) * mtrxr .* param.coil_r(1, k)) ./ tmp1;
        [K, E] = ellipke(kk);
        psi = psi + 2 * pi * Mu * param.coil_Ic(1, k) ./ (pi * (kk).^0.5) .* sqrt(param.coil_r(1, k) .* mtrxr) .* ...
            ((1 - kk / 2) .* K - E);
    end

end
