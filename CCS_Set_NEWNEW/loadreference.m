
%% loadreference!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadreference!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadreference!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadreference!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function REF = loadreference(PARAM)

    REF.Flux = dlmread([PARAM.input_file_directory '/FluxProfile_2D.txt']);

    if PARAM.IUTST==5
        REF.Flux =REF.Flux-0.0042459;
    end
    
    REF.Flux = REF.Flux / 2./pi;
    [znum, rnum] = size(REF.Flux);

    %% Here, the geometory should be derived from input file
    zmin = -0.9985;
    zmax = 0.9985;
    rmin = 0.108150;
    rmax = 0.6940;
    delz = (zmax-zmin)/(znum-1);
    delr = (rmax-rmin)/(rnum-1);

    for i=1:znum
        REF.Z(i) = zmin + delz * (i-1);
    end
    for i=1:rnum
        REF.R(i) = rmin + delr * (i-1);
    end
    REF.Flux(abs(REF.Z)>0.285, REF.R > 0.5985) = 0;
end
%% loadreference kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
