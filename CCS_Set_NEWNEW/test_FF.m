FF = load("FF.mat").FF;
FF
FF(end + 1:end + sum(9)) = 0.0;
FF

if 0
    FF(end + 1) = -66410 * 4.0 * pi * 1.0e-7;
end

FF
