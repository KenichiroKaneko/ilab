function dispSensorPosition(PARAM, SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, REF)
    figure()
    hold on
    axis equal
    VV = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);
    axis([0 1 -1.2 1.2])
    plot(VV(:, 1), VV(:, 2), '-k');
    xlabel("R[m]");
    ylabel("Z[m]");
    title("Magnetic Probe");

    if SENSOR_TPRB.NUM > 0
        scatter(SENSOR_TPRB.R, SENSOR_TPRB.Z, "r", "filled");
    end

    if SENSOR_NPRB.NUM > 0
        scatter(SENSOR_NPRB.R, SENSOR_NPRB.Z, "r");
    end

    figure()
    hold on
    axis equal
    VV = dlmread([PARAM.temporary_file_directory '/VacuumVesselMeshPoints.txt']);
    axis([0 1 -1.2 1.2])
    plot(VV(:, 1), VV(:, 2), '-k');
    xlabel("R[m]");
    ylabel("Z[m]");
    title("Flux Loop");
    scatter(SENSOR_FLXLP.R, SENSOR_FLXLP.Z, "b");

    for i = 1:PARAM.CCS
        % scatter(CCSDAT.RCCN(i, :), CCSDAT.ZCCN(i, :), 'mo', "filled");
    end

    v = linspace(-20, 20, 11);

    if REF ~= 0
        contour(REF.R, REF.Z, REF.Flux, '--k', "LineWidth", 1); % ????
    end

    % error('error description', A1)

    % if SENSOR_NPRB.NUM > 0
    %     legend("B-tangencial", "B-norm", "Flux", "CCSNode", "Ref-flux");
    % else
    %     legend("B", "Flux", "CCSNode", "Ref-flux");
    % end

    % title("Sensor positions and Reference");
end
