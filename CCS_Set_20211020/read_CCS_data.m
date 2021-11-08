% need shot_no
function [time, Plasma_Current, Coil_Current, Flux_Loop, Magnetic_Probe, Magnetic_Probe_Lowpass] = read_CCS_data(foldername, shotnum)

    foldername = '180515';
    shotnum = '010';

    lowpass_freq = 1000 * 1000;

    % close all

    %NI Jackie NJ lizzie ICS Will??Ahdr???w?b?_?[?Adat??dat
    %?S???????????????A???????R?????g????????L????????????????B

    %?S??I??ICS???????????X????K?v??????B??????X?_??H
    % U:~????????????A?l?b?g???[?N?h???C?u?o?R?????????
    % ??????g??????????AU:??UTST?h???C?u?????дн???Afoldrrname??shotno?????дн???

    %U:~???O???????s??A??????h?L???????g?i?R?[?h???????????????t?H???_?j??_?E?????[?h??????????????????B

    hdr_NJ = [foldername, '/NJ', shotnum, '.hdr'];
    dat_NJ = [foldername, '/NJ', shotnum, '.dat'];

    % ######################
    % ## READ HEADER DATA ##
    % ######################

    % read first line
    temp_NI = dlmread(hdr_NJ, '\t', [0 0 0 1]);
    ver_NI = temp_NI(1, 1);
    nch_NI = temp_NI(1, 2);
    % % READ headr
    fileID = fopen(hdr_NJ);
    header_NI = textscan(fileID, '%u %s %u %f64 %f64 %f64 %f64', 'HeaderLines', 1);
    fclose(fileID);
    chno = vertcat(header_NI{1});
    name = vertcat(header_NI{2});
    n = vertcat(header_NI{3});
    dt = vertcat(header_NI{4});
    range = vertcat(header_NI{6});
    gain = vertcat(header_NI{7});

    % pause

    % list_size = size(name);
    % pause
    % for i = 1 : list_size
    %     fprintf('%d %s %d \n', i, cell2mat(name(i)), gain(i))
    % end

    % ##########################
    % ###  READ BINARY DATA  ###
    % ##########################

    n = n(1);
    fid_NI = fopen(dat_NJ, 'r', 'b');
    data = horzcat(fread(fid_NI, [n, nch_NI], 'int16', 'b'));
    fclose(fid_NI);

    % #############################################
    % ###   CALIBLATE THE BINARY DATA           ###
    % ###  BY USING INFORMATION OF HEADER FILE  ###
    % #############################################

    time = zeros(1, 30000);
    dt = dt(1);
    I = 1:n;
    time(I) = double(I) * dt / double(1000);
    data_binary = zeros(n, nch_NI);

    for I = 1:nch_NI

        for J = 1:n
            data_binary(J, I) = range(I) / 32768 * gain(I) * data(J, I);
        end

    end

    % M?~N?z??@M?F?`?????l?????AN?F?v???_?i15ms?A30000?_?j
    % ##### input plsma current current data #####
    Plasma_Current = zeros(1, 30000);
    Plasma_Current(1, :) = data_binary(:, 1);
    % figure('Name','Plasma Current','NumberTitle','off')
    % plot(time, Plasma_Current)

    % ##### input coil current data #####
    Coil_Current = zeros(5, 30000);
    % figure('Name','Coil Current','NumberTitle','off')
    for i = 2:9

        if i == 9 % 16????PF4L_feed???i?[????????
            Coil_Current(8, :) = data_binary(:, 16);
            %         subplot(8, 1, 8); plot(time, Coil_Current(8, :))
            %         legend(name(16));
        else
            Coil_Current(i - 1, :) = data_binary(:, i);
            %         subplot(8, 1, i - 1); plot(time, Coil_Current(i - 1, :))
            %         legend(name(i));
        end

    end

    % ##### input flux loop data #####
    Flux_Loop = zeros(35, 30000);
    % figure('Name','Flux Loop Signal (Gain multiplied (here, 1 or -1))','NumberTitle','off')

    for i = 18:53

        if and(19 <= i, i <= 22)
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif and(27 <= i, i <= 29)
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif and(31 <= i, i <= 33)
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif and(37 <= i, i <= 38)
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif and(41 <= i, i <= 42)
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif and(44 <= i, i <= 45)
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif i == 48
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif and(50 <= i, i <= 51)
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif i == 53
            Flux_Loop(i - 18, :) = (-1) * data_binary(:, i);
        elseif i == 18 % 18????FL25???i?[????????
            Flux_Loop(25, :) = data_binary(:, i);
            %         subplot(5, 7, 25); plot(time, Flux_Loop(25, :))
            %         legend(name(i));
        elseif i == 43
            continue
        else
            Flux_Loop(i - 18, :) = data_binary(:, i);
            %         subplot(5, 7, i - 18); plot(time, Flux_Loop(i - 18, :))
            %         legend(name(i));
        end

    end

    % ##### input magnetic probe data #####
    Magnetic_Probe = zeros(40, 30000);
    Magnetic_Probe_FT = zeros(40, 30000);
    Magnetic_Probe_FT_LP = zeros(40, 30000);
    Magnetic_Probe_Lowpass = zeros(40, 30000);
    freq = fftshift((0:length(Magnetic_Probe(1, :)) -1) * 2000000 / length(Magnetic_Probe(1, :))) - 1000000;
    % figure('Name','Magnetic Probe Signal (Gain multiplied (here, 1 or -1))','NumberTitle','off')

    for i = 54:93

        if i == 74
            Magnetic_Probe(i - 53, :) = -data_binary(:, i);
        elseif and(79 <= i, i <= 82)
            Magnetic_Probe(i - 53, :) = -data_binary(:, i);
        elseif and(84 <= i, i <= 86)
            Magnetic_Probe(i - 53, :) = -data_binary(:, i);
        elseif i == 92
            Magnetic_Probe(i - 53, :) = -data_binary(:, i);
        else
            Magnetic_Probe(i - 53, :) = data_binary(:, i);
        end

    end

    for i = 54:93
        %     Magnetic_Probe(i - 53, :) = data_binary(:, i);
        Magnetic_Probe_FT(i - 53, :) = fft(Magnetic_Probe(i - 53, :));
        Magnetic_Probe_FT_lp = Magnetic_Probe_FT(i - 53, :);
        Magnetic_Probe_FT_lp(freq > lowpass_freq) = 0;
        Magnetic_Probe_FT_lp(freq < -lowpass_freq) = 0;
        Magnetic_Probe_FT_LP(i - 53, :) = Magnetic_Probe_FT_lp;
        Magnetic_Probe_Lowpass(i - 53, :) = ifft(Magnetic_Probe_FT_LP(i - 53, :));
        %     subplot(5, 8, i - 53); plot(time, Magnetic_Probe_Lowpass(i - 53, :))
        %     legend(name(i));
    end

end
