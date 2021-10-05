function data = UTSTdataload(dir, date, shotnum, type)
    p = join([dir 'utst/' num2str(date, '%06d') '/' type num2str(shotnum, '%03d') '.hdr'], "");

    if (exist(p, 'file') == 2)
        file = fopen(p, 'r');
    else
        disp(p);
        data = [];
        return;
    end

    version = fscanf(file, '%d', [1, 1]); % % 1
    num = fscanf(file, '%d', [1, 1]); % % ninum
    data_s = struct('ch', cell(1), 'name', cell(1), 'length', cell(1), 'freq', cell(1), 'offset', cell(1), 'range', cell(1), 'factor', cell(1), ...
        'raw', cell(1), 'sig', cell(1), 'time', cell(1));
    data = repmat(data_s, num, 1);

    for i = 1:num
        data(i).ch = fscanf(file, '%d', [1, 1]); % % ch_number
        data(i).name = fscanf(file, '%s', [1, 1]); % % ch_name
        data(i).length = fscanf(file, '%d', [1, 1]); % % ch_length
        data(i).freq = fscanf(file, '%f', [1, 1]) * 1e-6; % % ch_freq_in_microsec
        data(i).offset = fscanf(file, '%f', [1, 1]); % % ch_offset
        data(i).range = fscanf(file, '%f', [1, 1]); % % ch_fullrange
        data(i).factor = fscanf(file, '%f', [1, 1]); % % ch_factor_to_be_multiplied
    end

    fclose(file);

    p = join([dir 'utst/' num2str(date, '%06d') '/' type num2str(shotnum, '%03d') '.dat'], "");

    switch version
        case 1

            if (exist(p, 'file') == 2)
                file = fopen(p, 'r');
            else
                disp(p);
                data = [];
                return;
            end

        case 2

            if (exist(p, 'file') == 2)
                file = fopen(p, 'r', 'b');
            else
                %disp(['File not found! (' num2str(date,'%06d') '/' num2str(shotnum,'%03d') '.dat)']);
                data = [];
                return;
            end

    end

    %if file==-1
    %    disp(['File not found! (' num2str(date,'%06d') '/' num2str(shotnum,'%03d') '.dat)']);
    %    data=[];
    %    return;
    %else
    for i = 1:num
        data(i).raw = fread(file, [data(i).length, 1], 'int16');
        data(i).sig = data(i).raw / 2^15 * data(i).range * data(i).factor;
        data(i).time = 0:data(i).freq:data(i).freq * (data(i).length - 1);
    end

    fclose(file);
    % %end
