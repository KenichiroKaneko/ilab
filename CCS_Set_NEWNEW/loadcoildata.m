%% loadcoildata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadcoildata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadcoildata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadcoildata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function ExtCOIL = loadcoildata(PARAM)

    %% Here, the geometory should be derived from input file
    ExtCOIL.NUM = 10;
    ExtCOIL.NAME = ['EFL', 'EFU', 'PF1L', 'PF1U', 'PF2L', 'PF2U',  'PF3L', 'PF3U', 'PF4L', 'PF4U'];
    ExtCOIL.R = [0.80, 0.80, 0.20, 0.20, 0.665, 0.665, 0.750, 0.750, 0.685, 0.685];
    ExtCOIL.Z = [-1.07, 1.07, -1.10, 1.10, -0.80, 0.80, -0.675, 0.675, -0.50, 0.50];
    ExtCOIL.C = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    ExtCOIL.N = [200, 200, 8, 8, 3, 3, 8, 8, 3, 3]; 
    ExtCOIL.I = [0.28, 0.28, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
end
%% loadcoildata kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
