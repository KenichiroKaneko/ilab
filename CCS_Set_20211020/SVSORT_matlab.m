%function [SVS] = SVSORT_matlab(PARAM,N,SVO,NMX)
function [SVS] = SVSORT_matlab(PARAM, SVO)
    N = length(SVO);
    SVS = zeros(N);
    fid13 = fopen([PARAM.temporary_file_directory '/@SV_SORT.txt'], 'w'); %13
    fid14 = fopen([PARAM.temporary_file_directory '/@SV_STEEP.txt'], 'w'); %14
    SVMX = -1.0D10;

    % SVSの最大値を求めている
    for K = 1:N
        SVS(K) = SVO(K);

        if (SVS(K) > SVMX)
            SVMX = SVS(K);
        else
        end

    end

    % fprintf('%s %d\r\n', 'SVMX =', SVMX);

    % SVSを大きい順に並べ替えている
    for J = 2:N
        A = SVS(J);

        for I = J - 1:-1:1

            if (SVS(I) >= A)
                break
            end

            SVS(I + 1) = SVS(I);
            I = 0;
        end

        SVS(I + 1) = A;
    end

    % SVSをSVSの最大値で割って規格化している
    SVS(1:N) = SVS(1:N) / SVMX;

    for K = 1:N
        fprintf(fid13, '%d %d\r\n', K, SVS(K));
    end

    % SVSの変化の割合をK,STEEPに保存している
    for K = 1:N - 1
        EE = SVS(K + 1) / SVS(K);
        STEEP = log(EE);
        fprintf(fid14, '%d %d\r\n', K, STEEP);
    end

    % save("vars_SVSORT");

end
