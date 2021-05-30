function [SENSOR_FLXLP, SENSOR_TPRB] = loadRealCCSdata()

    foldername = '180515';
    shotnum = '010';
    shutter = 17700;
    [time, Plasma_Current, Coil_Current, Flux_Loop, Magnetic_Probe, Magnetic_Probe_Lowpass] = read_CCS_data(foldername, shotnum);
    size(time)
    size(Plasma_Current)
    size(Coil_Current)
    size(Flux_Loop)
    size(Magnetic_Probe)
    size(Magnetic_Probe_Lowpass)

    SENSOR_FLXLP = Flux_Loop(shutter);
    SENSOR_TPRB = Magnetic_Probe(shutter);

    % 0.0005~15msecまで 0.0005msecごとにある
    figure()
    hold on
    plot(time, Plasma_Current)
    scatter(time(shutter), Plasma_Current(shutter))

end
