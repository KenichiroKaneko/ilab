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

fid61 = fopen('CCS_temporary/Parameters_FL_clockwise_180515010_t8000_EFkaizen.txt', 'r');
fid68 = fopen('CCS_temporary/Parameters_MP_clockwise_180515010_t8000_EFkaizen.txt', 'r');
fid62 = fopen('CCS_temporary/@UTST_SenPos.txt', 'w');
IMAX_b = 39;
IMAX_f = 33;

for I = 1:IMAX_f
    temp = textscan(fid61, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R = str2double(temp_m(1));
    Z = str2double(temp_m(2));
    fprintf(fid62, '%d %d\r\n', R, Z);
    fprintf(fid62, '%d %d\r\n', R, -Z);
end

for I = 1:IMAX_b
    temp = textscan(fid68, '%s', 1, 'Delimiter', '\n');
    temp_m = strsplit(temp{1}{1});
    R = str2double(temp_m(1));
    Z = str2double(temp_m(2));
    %PSI = str2double(temp_m(3));
    %     BZ = str2double(temp_m(4));
    BR = str2double(temp_m(5)); %      READ(61,*),R,Z,PSI,BZ,BR

    BZ = BZ_Distribution(1, I);

    %
    E = abs((R - RSENS) / RSENS);
    LG = LG + 1;
    LLG = LG;

    if (LG == LGB) % 1?��?��?��?��?���??��?��Z?��b?��g?��?��?��?��
        LG = 0;
    end

    MG = MG + 1;
    MMG = MG;

    if (MG == MGB)
        MG = 0;
    end

    %%C      IF(LLG.NE.LGBH) GOTO 200
    if and (E < 1.0D - 5, MMG ~= MGBH)
        continue
    elseif and (E >= 1.0D - 5, LLG ~= LGBH)
        continue

    elseif and(LWALL > 0, and (R < 0.12, abs(Z) < 0.9)); % ?��ǉ�?��?��?��̃Z?��?��?��T?��[?��?��?��?��2016.9.5ushiki
        continue
    end

    RR = R;
    ZZ = -Z;
    REND = RR - BR * FACTR;
    ZEND = ZZ + BZ * FACTR;

    II = II + 1;

    RS(NFLX + 2 * II - 1) = R;
    RS(NFLX + 2 * II) = R;
    ZS(NFLX + 2 * II - 1) = Z;
    ZS(NFLX + 2 * II) = -Z;
    Z1 = Z; BR1 = BR;
    Z2 = -Z; BR2 = -BR;
    TET1 = atan2(BZ, BR1);
    TET2 = atan2(BZ, BR2);
    % ?��?��?��͐�?��̐ڐ�?��?��?��?��?��̎�?��?��?��?��E?��?��?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I?��I
    TET(NFLX + 2 * II - 1) = TET1; %  !
    TET(NFLX + 2 * II) = TET2; %  !?��@
    ITYPE(NFLX + 2 * II - 1) = 1;
    ITYPE(NFLX + 2 * II) = 1;
    BB1 = BR1 * cos(TET1) + BZ * sin(TET1);
    BB2 = BR2 * cos(TET2) + BZ * sin(TET2);
    XBR1 = BB1 * cos(TET1);
    XBZ1 = BB1 * sin(TET1);
    XBR2 = BB2 * cos(TET2);
    XBZ2 = BB2 * sin(TET2);

    EPS1R = abs((XBR1 - BR1) / BR1);
    EPS1Z = abs((XBZ1 - BZ) / BZ);

    EPS2R = abs((XBR2 - BR2) / BR2);
    EPS2Z = abs((XBZ2 - BZ) / BZ);

    TPRB(2 * II - 1) = BB1;
    TPRB(2 * II) = BB2;

end %200
