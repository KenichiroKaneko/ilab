% vars = load('vars_loadLizzieData');
% Coil_Current = vars.Coil_Current;
% time = vars.time;
% time_CCS = vars.time_CCS;

% figure
% plot(time, smoothdata(Coil_Current(2 * i - 1, :), 'gaussian', 200))

% f = smoothdata(Coil_Current(2 * i - 1, :), 'gaussian', 200);
% f(time_CCS / (0.5 * 0.001))

for i = 3:10

    if rem(i, 2)
        i = i + 1;
    else
        i = i - 1;
    end

    i
    % ExtCOIL.I(i)
end

% 3 4 5 6 7 8 9 10
% 4 3 6 5 8 7 10 9
for i = 1:8

    if rem(i, 2)
        i = i + 3;
    else
        i = i + 1;
    end

    i
    % ExtCOIL.I(i)
end
