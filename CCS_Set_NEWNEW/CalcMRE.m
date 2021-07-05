function CalcMRE(PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, A, ss, vv, uu, X, FC, FF)
    % 各センサーの相対誤差を計算

    NFLX = SENSOR_FLXLP.NUM;
    NAPB = SENSOR_TPRB.NUM + SENSOR_NPRB.NUM;
    NCCN = sum(CCSDAT.NCCN);

    NFLXl = 1;
    NFLXr = NFLX;
    NAPBl = 1 + NFLXr;
    NAPBr = NFLXr + NAPB;
    NCCSl = 1 + NAPBr;
    NCCSr = NAPBr + NCCN;

    % RE = ((A * X' - FC') ./ FC').^2;
    RE = abs(A * X' - FF');

    MREFlux = norm(RE(NFLXl:NFLXr)) / NFLX;
    MREAPB = norm(RE(NAPBl:NAPBr)) / NAPB;
    MRECCS = norm(RE(NCCSl:NCCSr)) / NCCN;

    MREs = [MREFlux, MREAPB, MRECCS];
    % disp(MREs);
    txt1 = sprintf("Flux Sensor   Num:%4d  MRE:%.3e\n", NFLX, MREFlux);
    txt2 = sprintf("Probe Sensor  Num:%4d  MRE:%.3e\n", NAPB, MREAPB);
    txt3 = sprintf("CCS Node      Num:%4d  MRE:%.3e\n", NCCN, MRECCS);
    disp(txt1 + txt2 + txt3)

    % ふし目を探索
    R = SENSOR_FLXLP.R;
    Z = SENSOR_FLXLP.Z;
    count_FLUX = 0;

    for i = 1:SENSOR_FLXLP.NUM - 1

        if (R(i) ~= R(i + 1) && Z(i) ~= Z(i + 1)) || (i == 1) || (i == NAPB)
            count_FLUX = count_FLUX + 1;
            corner_FLUX(count_FLUX) = i;
        end

    end

    R = SENSOR_NPRB.R;
    Z = SENSOR_NPRB.Z;
    count_NPRB = 0;

    for i = 1:SENSOR_NPRB.NUM - 1

        if (R(i) ~= R(i + 1) && Z(i) ~= Z(i + 1)) || (i == 1) || (i == SENSOR_NPRB.NUM - 1)
            count_NPRB = count_NPRB + 1;
            corner_NPRB(count_NPRB) = NFLX + i
        end

    end

    if CONFIG.ShowFig
        figure()
        hold on
        plot(NFLXl:NCCSr, A * X');
        plot(NFLXl:NCCSr, FF);
        % plot(NFLXl:NCCSr, RE);

        x = [NFLXr, NFLXr, NAPBl, NAPBr, NCCSl, NCCSr];
        scatter(x, FF(x), "m*")

        % for i = 1:count_FLUX
        %     scatter(corner_FLUX(i), FC(corner_FLUX(i)), 'm*')
        % end

        % for i = 1:count_NPRB
        %     scatter(corner_NPRB(i), FC(corner_NPRB(i)), 'm*')
        % end

        title('各センサー位置の再構成結果');
        legend('A*X', 'FF');
    end

    % save("vars_inMRE  ");
    % error('f', A1)
end
