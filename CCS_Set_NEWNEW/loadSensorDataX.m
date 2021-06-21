function [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R] = loadsensordata(PARAM, CONFIG)

    CCS_R = 0;
    CCS_Z = 0;

    % データの読み込み
    % 0:.mat or 1:.txt
    if CONFIG.FileExt
        % txtファイルでの読み込み

        if CONFIG.DataType
            % 0:GS or 1:exp
            % 実際のセンサー位置で再構成するコード、入力は末尾がRと、実験データが対応(例)UTST_numel_2033R、180515_010_t9650
            [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R] = loadRealsensordata(PARAM);

            % 法線方向、接線方向を考慮して読み込むためのコード、入力は末尾がRと、実験データが対応
            % [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z,CCS_R] = loadRealsensordataTN(PARAM);
        else
            % オリジナルのコード、入力は、末尾にRがないもの、UTST_numel_5,9が対応
            [SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCS_Z, CCS_R] = loadsensordata(PARAM);
        end

    else
        % matファイルでの読み込み
        vars = load(PARAM.input_file_directory + "/merged.mat");
        R = vars.r_CCS;
        Z = vars.z_CCS;
        BR = vars.Br_CCS;
        BZ = vars.Bz_CCS;
        PSI = vars.psi_CCS;

        % センサー位置でデータを間引く
        len = length(R);
        index = [1:1:len];
        R = R(index);
        Z = Z(index);
        BR = BR(index);
        BZ = BZ(index);
        PSI = PSI(index);
        len = length(R);

        SENSOR_TPRB.R = R;
        SENSOR_TPRB.Z = Z;
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
        RR = R;
        ZZ = Z;
        B = sqrt(BR.^2 + BZ.^2);
        SENSOR_TPRB.TNPRB = sqrt(BR.^2 + BZ.^2);
        SENSOR_TPRB.NUM = len;
        disp(['Number of TPRB =  ' num2str(SENSOR_TPRB.NUM)]);

        SENSOR_FLXLP.R = R;
        SENSOR_FLXLP.Z = Z;

        if CONFIG.DevideFlux
            SENSOR_FLXLP.FLXLP = PSI ./ (pi .* R .* R);
        else
            SENSOR_FLXLP.FLXLP = PSI;
        end

        SENSOR_FLXLP.TET = zeros(1, len);
        SENSOR_FLXLP.ITYPE = zeros(1, len);
        SENSOR_FLXLP.NUM = len;
        disp(['Number of FLXLP = ' num2str(SENSOR_FLXLP.NUM)]);

        save("vars_inloadX")
    end

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

    end

    % CCS面の自動決定R方向
    if CONFIG.DetermineCCSRPos
        CCS_R = CalcPlasmaCenter(PARAM, CONFIG, CCS_Z);
    end

end
