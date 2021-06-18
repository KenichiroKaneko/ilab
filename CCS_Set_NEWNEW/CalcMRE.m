function CalcMRE(PARAM, CONFIG, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, A, ss, vv, uu, X, FC)
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

    RE = (A * X' - FC) ./ FC;

    MREFlux = norm(RE(NFLXl:NFLXr)) / NFLX;
    MREAPB = norm(RE(NAPBl:NAPBr)) / NAPB;
    MRECCS = norm(RE(NCCSl:NCCSr)) / NCCN;

    MREs = [MREFlux, MREAPB, MRECCS];
    % disp(MREs);
    txt1 = sprintf("Flux Sensor   Num:%4d  MRE:%.3e\n", NFLX, MREFlux);
    txt2 = sprintf("Probe Sensor  Num:%4d  MRE:%.3e\n", NAPB, MREAPB);
    txt3 = sprintf("CCS Node      Num:%4d  MRE:%.3e\n", NCCN, MRECCS);
    disp(txt1 + txt2 + txt3)

    save("vars_inMRE");
    % error('f', A1)
end
