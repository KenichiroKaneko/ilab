function dispSensorPosition(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, REF)
    figure()
    hold on

    if SENSOR_TPRB.NUM > 0
        scatter(SENSOR_TPRB.R, SENSOR_TPRB.Z, "r", "filled");
    end

    if SENSOR_NPRB.NUM > 0
        scatter(SENSOR_NPRB.R, SENSOR_NPRB.Z, "r");
    end

    scatter(SENSOR_FLXLP.R, SENSOR_FLXLP.Z, "filled");

    for i = 1:PARAM.CCS
        scatter(CCSDAT.RCCN(i, :), CCSDAT.ZCCN(i, :), 'mo', "filled");
    end

    v = linspace(-20, 20, 11);
    contour(REF.R, REF.Z, REF.Flux, '--k', "LineWidth", 1); % ????
    xlabel("R[m]");
    ylabel("Z[m]");
    hold off
    legend("B-tangencial", "B-norm", "Flux", "CCSNode");
    title("Sensor positions and Reference");
end
