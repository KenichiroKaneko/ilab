%% loadinputfile!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadinputfile!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadinputfile!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% loadinputfile!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function PARAM = loadinputfile(filename)

    %inputdata0=fileread('CCS_input\input.txt');
    inputdata0=fileread(['CCS_input\' filename]);
    inputdata=strsplit(inputdata0,{'\n','\t','\r'});

    % file direcory
    PARAMNUM=1;
    PARAM.input_file_directory = inputdata{PARAMNUM}; PARAMNUM=PARAMNUM+2;
    
    PARAM.temporary_file_directory = inputdata{3}; PARAMNUM=PARAMNUM+2;
    
    PARAM.output_file_directory = inputdata{5}; PARAMNUM=PARAMNUM+2;

    % CCS parameters
    PARAM.IUTST =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.GETA_YN = str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.NONC =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.SIGM =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.SEED =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.CCS =     str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.IDECCS =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;

    for i=1:PARAM.CCS
        PARAM.R0(i) =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.Z0(i)=    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.RR(i) =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.CAPPER(i) =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.NE(i) =   str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.MSEC(i,1) =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.MSEC(i,2) =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.MSEC(i,3) =    str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
        PARAM.SOU(i) =      str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    end
    
    PARAM.IPCONST =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.LCOND =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.ITSKP =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.ITRNC =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.IDCN =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.KUP0 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.MTS =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.LXMAP =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.PLOTW =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.SERX =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.LOCAT =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.AAA1 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.BBB1 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.AAA2 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.BBB2 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.statR =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.endR =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.I0 =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
    PARAM.TRUTOT =  str2double(inputdata{PARAMNUM}); PARAMNUM=PARAMNUM+2;
end
%% loadinputfile kokomade!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
