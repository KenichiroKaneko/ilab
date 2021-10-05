dir = "./";
date = 210906;
shotnum = 7;
s = 9510;

data_NJ = UTSTdataload(dir, date, shotnum, 'NJ');

% FLのプロット
f4 = figure('Name', 'FL');

FL = [];
MP = [];

for i = 1:35
    ch = i + 18;

    if (ch == 25)
        ch = 18;
    end

    data_NJ(ch).time(end - 24:end) = [];
    data_NJ(ch).raw(end - 24:end) = [];
    data_NJ(ch).sig(end - 24:end) = [];
    offset = mean(data_NJ(ch).sig(50:round(450 / data_NJ(ch).freq / 1e6)));
    % 積分
    % data_NJ(ch).sigint = cumtrapz(data_NJ(ch).sig - offset) * data_NJ(ch).freq;

    figure(f4);
    hold on
    % M = movmean(data_NJ(ch).sig, 10); % 移動平均
    % l = 1; r = 200;
    % M = M - sum(data_NJ(ch).sig(l:r)) / (r - l + 1); % オフセット調整
    % plot(data_NJ(ch).time, M)
    FL(i) = data_NJ(ch).sig(time2Frame(s));
end

for i = 1:40
    ch = i + 53;
    M = movmean(data_NJ(ch).sig, 10); % 移動平均
    l = 1; r = 200;
    data_NJ(ch).sig = data_NJ(ch).sig - sum(data_NJ(ch).sig(l:r)) / (r - l + 1); % オフセット調整
    data_NJ(i).sigint = cumtrapz(data_NJ(i).sig - offset) * data_NJ(i).freq;
    % plot(data_NJ(ch).time, M)
    MP(i) = data_NJ(ch).sig(time2Frame(s));
end

out = [MP FL];
figure()
plot(1:length(out), out);

% 15000µsecを30000frameで観測している
function frame = time2Frame(time)
    frame = time * 2;
    return;
end
