tic

% 諸変数の定義
EPS = 1.0e-8; % 許容誤差
Mu = 4 * pi * 1.0e-7; % 真空の透磁率
kappa = 1.38e-23; % ボルツマン定数

% Sw_initを受け取る
fprintf("usage : gs [switch of initial condition]\n");
fprintf(" switch == 1  ... calculate external flux\n");
fprintf(" switch == 2  ... use previous external flux\n");
fprintf(" switch == 3  ... use previous equilibrium\n");
prompt = "input switch number";
% Sw_init = input(prompt);
% global Sw_init;
param.Sw_init = 1;

% ファイル名や定数の定義
% datafile_dir = "datafiles_CCS03c_finemesh/";
% output_dir = "outputfiles_03c/";
% datafile_dir = "datafiles_CCS04c_finemesh/";
% output_dir = "outputfiles_04c/";
% datafile_dir = "datafiles_CCS05c_double_finemesh/";
% output_dir = "outputfiles_05c/";
% datafile_dir = "datafiles_test/";%入力データ
% output_dir = "outputfiles_test/";% 出力データフォルダ

% 大量の数値解を作成するためのファイル名一覧
% filenames = ["z2033_r602" "z1000_r602" "z0960_r602" "z0920_r602" ...
    %             "z0880_r602""z0840_r602""z0800_r602""z0760_r602" ...
%             "z0720_r602""z0680_r602""z0640_r602""z0600_r602" ...
%             "z0560_r602""z0520_r602""z0480_r602""z0440_r602" ...
%             "z0400_r602""z0360_r602""z0320_r602""z0280_r602" ...
%             "z0240_r602""z0200_r602""z0160_r602""z0120_r602" ...
%             "z0080_r602"];

% filenames = ["z0240_r602" "z0200_r602" "z0160_r602" "z0120_r602" "z0080_r602"];
filenames = ["z2033_r602"];

for fileNum = 1:length(filenames)
    disp(["****  " + filenames(fileNum) + "  " + fileNum + "/" + length(filenames) + "  ****"])

    datafile_dir = "INPUT/" + filenames(fileNum) + "/";
    output_dir = "OUTPUT/" + filenames(fileNum) + "/";

    param.file0 = datafile_dir + "def_mesh.dat";
    param.file1 = datafile_dir + "def_slow.dat";
    param.file2 = datafile_dir + "def_gs.dat";
    param.file3 = datafile_dir + "def_fast.dat";

    % ファイルから値の読み込み
    A = load(param.file0);
    param.Nz = A(1, 1);
    param.Nr = A(2, 1);
    % param.Nz = 2033;
    % param.Nr = 602;
    param.Nz_orig = param.Nz;
    param.Nr_orig = param.Nr;
    param.zmin = A(3, 1);
    param.zmax = A(4, 1);
    param.rmin = A(5, 1);
    param.rmax = A(6, 1);
    param.z1 = A(7, 1);
    param.z2 = A(8, 1);
    param.z3 = A(9, 1);
    param.z4 = A(10, 1);
    param.rmirror = A(11, 1);
    % param.zlimiter = A(12, 1);
    % param.rlimiter = A(13, 1);
    param.routine_num = 0;
    param.zlimiter = 0;
    param.rlimiter = 0.1656627;

    % ファイルから値の読み込み
    B = load(param.file2);
    param.alpha = B(1, 1);
    param.beta = B(2, 1);
    param.P0 = B(3, 1); % 0~500~
    param.Iplasma = B(4, 1); % 大きさ
    param.Itfc = B(5, 1); % あまり変わらない
    param.errormax = B(6, 1);
    param.itermax = B(7, 1);
    param.omega = B(8, 1);
    param.RBt = Mu * param.Itfc / 2 / pi;

    % 初期状態を設定
    param = define_params(param, param.Nz, param.Nr);

    %% DEFINE OMEGA
    param.omega = 1.945;

    %% 荒いメッシュで計算
    [param, psi, mu_jt, p, Ip] = sub_routine_init(param, 4, 4);
    [param, psi, mu_jt, p, Ip] = sub_routine(param, 4, 4, psi, mu_jt, p, Ip);
    [param, psi, mu_jt, p, Ip] = sub_routine(param, 2, 2, psi, mu_jt, p, Ip);
    [param, psi, mu_jt, p, Ip] = sub_routine(param, 1.5, 1.5, psi, mu_jt, p, Ip);

    % [param, psi, mu_jt, p, Ip] = sub_routine_init(param, 1, 1);
    % [param, psi, mu_jt, p, Ip] = sub_routine(param, 1, 1, psi, mu_jt, p, Ip);

    %% メインのルーチン用に変形
    param.Nz = param.Nz_orig;
    param.Nr = param.Nr_orig;
    psi = imresize(psi, [param.Nz, param.Nr]);
    mu_jt = zeros(param.Nz, param.Nr);
    p = zeros(param.Nz, param.Nr);
    Ip = zeros(param.Nz, param.Nr);
    param = define_params(param, param.Nz, param.Nr);
    psi = psi_round_revive(psi, param);

    %% メインのルーチン
    for i = 1:10000000
        [psi_bar, mu_jt, p, Ip] = cal_jt(param, psi, mu_jt, p, Ip);
        % 真空容器の中にコイルがある時
        % mu_jt = addPFcurrent(param,psi,mu_jt);
        format long
        [error, psi] = cal_flux_C(param, psi, mu_jt);

        if (rem(i, 100) == 0)
            disp([num2str(i) ':' num2str(error)])
        end

        if error < param.errormax
            param.routine_num = param.routine_num + i;
            disp(['Convergence! Iteration number is ' num2str(param.routine_num)]);
            break
        end

    end

    toc
    disp("fin")

    % WriteOut_data(psi, output_dir + "flux.txt");
    % WriteOut_data(mu_jt, output_dir + "mujt.txt");
    % WriteOut_data(Ip, output_dir + "rbt.txt");
    % WriteOut_data(p, output_dir + "p.txt");

    psi0 = zeros(param.Nz, param.Nr);
    [param, psi0] = cal_psi(param, psi0);

    % 変数の保存
    save(output_dir + "vars", "param", "psi", "mu_jt", "Ip", "p", "psi0");
    gs_graph(param, psi, p, Ip, mu_jt)
end

toc

%% psiを更新するルーチン。
function [error, psi] = cal_flux_C(param, psi, mu_jt)
    error = 0;
    vec3 = 1 / (param.delz^2);
    vec4 = 1 / (param.delz^2);

    for j = 2:param.Nr - 1
        radius = param.rmin + (j - 1) * param.delr;
        vec1 = 1 / (param.delr^2) * radius / (radius + 0.5 * param.delr);
        vec2 = 1 / (param.delr^2) * radius / (radius - 0.5 * param.delr);
        vec5 = vec1 + vec2 + vec3 + vec4;

        for i = 2:param.Nz - 1

            if param.psi_inside(i, j) == 1
                psi_old = psi(i, j);
                psi(i, j) = param.omega * (...
                    vec1 * psi(i, j + 1) + ...
                    vec2 * psi(i, j - 1) + ...
                    vec3 * psi(i + 1, j) + ...
                    vec4 * psi(i - 1, j) + ...
                    2 * pi * radius * mu_jt(i, j)) / vec5 + ...
                    (1 - param.omega) * psi_old;

                if (abs(psi(i, j) - psi_old) > error)
                    error = abs(psi(i, j) - psi_old);
                end

            end

        end

    end

end

%% ψの初期状態を置く
function [param, psi] = cal_psi(param, psi)
    Mu = 4 * pi * 1.0e-7; % 真空の透磁率
    mtrxz = (0:1:param.Nz - 1)' .* param.delz + param.zmin;
    mtrxr = (0:1:param.Nr - 1) .* param.delr + param.rmin;

    for k = 1:param.ncoil
        % (z*rの大きさの配列を作成)
        tmp1 = (param.coil_z(1, k) - mtrxz).^2 * ones(1, param.Nr) + ones(param.Nz, 1) * (param.coil_r(1, k) + mtrxr).^2;
        kk = (4 .* ones(param.Nz, 1) * mtrxr .* param.coil_r(1, k)) ./ tmp1;
        [K, E] = ellipke(kk);
        psi = psi + 2 * pi * Mu * param.coil_Ic(1, k) ./ (pi * (kk).^0.5) .* sqrt(param.coil_r(1, k) .* mtrxr) .* ...
            ((1 - kk / 2) .* K - E);
    end

end

%% psiの初期値を適当に決める
function [param, psi] = cal_initial_condition(param, psi)
    psi_axis = -1.75;
    psi_edge = 0;
    zc = (param.Nz - 1) / 2 * param.delz;
    rc = (param.Nr - 1) / 3 * param.delr;
    mtrxz = (0:1:param.Nz - 1)' .* param.delz;
    mtrxr = (0:1:param.Nr - 1) .* param.delr;
    radius = sqrt(((mtrxz - zc) ./ zc).^2 * ones(1, param.Nr) + ones(param.Nz, 1) * (mtrxr ./ rc).^2);
    psi_b = (psi_axis - psi_edge) * (1 - radius);
    % 2020/12/03
    psi_b(psi_b > 0) = 0;
    psi = psi - psi .* param.psi_inside + psi_b .* param.psi_inside;

end

%% 磁束の値から電流jtを求める
function [psi_bar, mu_jt, p, Ip] = cal_jt(param, psi, mu_jt, p, Ip)

    EPS = 1.0e-8; % 許容誤差
    Mu = 4 * pi * 1.0e-7; % 真空の透磁率

    psi_axis = min(psi(param.psi_inside == 1));
    psi_edge = min(psi(param.psi_round == 1));

    if (param.z_limiterpos * param.r_limiterpos ~= 0)
        psi_edge = min(psi_edge, psi(param.z_limiterpos, param.r_limiterpos));
    end

    % psi_edge = min(psi(66, :));

    if (psi_edge <= (psi_axis) + EPS)
        error("psi_edge <= psi_axis ...");
    end

    radius = param.rmin + ones(param.Nz, 1) * (0:1:param.Nr - 1) .* param.delr;
    psi_resion = param.psi_inside;
    psi_resion(psi > psi_edge | psi < psi_axis) = 0;
    psi_deff = psi_axis - psi_edge;
    ip = CALGP_t(param, (psi - psi_edge) .* psi_resion / psi_deff) / psi_deff;

    aaa = sum(((psi - psi_edge) .* psi_resion ./ psi_deff).^(2 * param.beta - 1) ./ radius ./ psi_deff, "all");
    bbb = sum(((psi - psi_edge) .* psi_resion ./ psi_deff).^(param.beta - 1) ./ radius ./ psi_deff, "all");
    ccc = sum(radius .* ip .* psi_resion, "all");

    aaa = aaa * (2 * pi * param.RBt^2 * param.beta * param.delr * param.delz);
    bbb = bbb * (-2 * pi * param.RBt^2 * param.beta * param.delr * param.delz);
    ccc = ccc * (2 * pi * Mu * param.delr * param.delz);
    ccc = ccc - param.Iplasma * Mu;

    if abs(aaa) < EPS
        param.gamma = -ccc / bbb;
    else
        abc = [aaa, bbb, ccc];
        sol = min(roots(abc));

        if isreal(sol)
            param.gamma = sol;
        else
            error("Solutions are complex.");
        end

    end

    radius = param.rmin + ones(param.Nz, 1) * (0:1:param.Nr - 1) .* param.delr;

    mu_jt = 2 * pi .* (radius .* Mu ...
        .* CALGP_t(param, (psi - psi_edge) ./ psi_deff) ./ psi_deff ...
        + CALIp_t(param, (psi - psi_edge) ./ psi_deff) ...
        .* CALGIp_t(param, (psi - psi_edge) ./ psi_deff) ./ psi_deff ./ radius);
    p = CALp_t(param, (psi - psi_edge) ./ psi_deff);
    Ip = CALIp_t(param, (psi - psi_edge) ./ psi_deff);
    mu_jt = mu_jt .* param.psi_inside;
    p = p .* param.psi_inside;
    Ip = Ip .* param.psi_inside;
    psi_bar = psi_axis / psi_edge;
end

function pp = CALp_t(param, psibar)
    pp = zeros(param.Nz, param.Nr);
    pp(psibar >= 0) = param.P0 .* power(psibar(psibar >= 0), param.alpha);
    pp(psibar < 0) = 0;
end

function gpp = CALGP_t(param, psibar)
    gpp = zeros(param.Nz, param.Nr);
    gpp(psibar >= 0) = param.alpha .* param.P0 .* power(psibar(psibar >= 0), param.alpha - 1);
    gpp(psibar < 0) = 0;
end

function ip = CALIp_t(param, psibar)
    ip = zeros(param.Nz, param.Nr);
    ip(psibar >= 0) = param.RBt .* (1 - param.gamma .* power(psibar(psibar >= 0), param.beta));
    ip(psibar < 0) = param.RBt;
end

function gip = CALGIp_t(param, psibar)
    gip = zeros(param.Nz, param.Nr);
    gip(psibar >= 0) = (-param.RBt * param.gamma * param.beta) .* power(psibar(psibar >= 0), param.beta - 1);
    gip(psibar < 0) = 0;
end

%% psiを更新するルーチン,配列を利用したものだがうまく動かない
function [error, psi] = cal_flux(param, psi, mu_jt)
    vec3 = 1 / (param.delz^2);
    vec4 = 1 / (param.delz^2);

    radius = param.rmin + ones(param.Nz, 1) * (0:1:param.Nr - 1) .* param.delr;
    vec1 = 1 / (param.delr^2) .* radius ./ (radius + 0.5 * param.delr);
    vec2 = 1 / (param.delr^2) .* radius ./ (radius - 0.5 * param.delr);

    vec5 = vec1 + vec2 + vec3 + vec4;

    psi_old = psi;
    psi_new = param.omega * (...
        vec1 .* circshift(psi, [0 -1]) ...
        +vec2 .* circshift(psi, [0 1]) ...
        +vec3 .* circshift(psi, [1 0]) ...
        +vec4 .* circshift(psi, [-1 0]) ...
        +2 * pi .* radius .* mu_jt) ./ vec5 ...
        +(1 - param.omega) * psi_old;

    psi = psi_old .* param.psi_round + psi_new .* param.psi_inside;

    C = max(abs(psi(param.psi_inside == 1) - psi_old(param.psi_inside == 1)));
    error = C;
end

%%
function mu_jt = addPFcurrent(param, psi, mu_jt)
    Mu = 4 * pi * 1.0e-7; % 真空の透磁率
    j = fix((param.coil_r(1, :) - param.rmin) ./ param.delr);
    r1 = param.coil_r(1, :) - (param.rmin + j * param.delr);
    i = fix((param.coil_z(1, :) - param.zmin) ./ param.delz);
    z1 = param.coil_z(1, :) - (param.zmin + i * param.delz);
    cur = Mu * param.coil_Ic(:)' ./ (param.delr^2 * param.delz^2);

    mu_jt(indexToMatrix(param, i + param.Nz .* j)) = mu_jt(indexToMatrix(param, i + param.Nz .* j)) + cur .* (param.delr - r1) .* (param.delz - z1);
    mu_jt(indexToMatrix(param, i + param.Nz .* (j + 1))) = mu_jt(indexToMatrix(param, i + param.Nz .* (j + 1))) + cur .* (r1) .* (param.delz - z1);
    mu_jt(indexToMatrix(param, i + 1 + param.Nz .* j)) = mu_jt(indexToMatrix(param, i + 1 + param.Nz .* j)) + cur .* (param.delr - r1) .* (z1);
    mu_jt(indexToMatrix(param, i + 1 + param.Nz .* (j + 1))) = mu_jt(indexToMatrix(param, i + 1 + param.Nz .* (j + 1))) + cur .* (r1) .* (z1);

    function [z, r] = indexToMatrix(param, num)
        z = fix(num ./ param.Nr);
        r = rem(num, param.Nr);

        if r == 0
            r = param.Nr;
        end

    end

end

% コマンドから実行する関数を選択するための関数
function [param, psi] = cal_initial_flux(param, psi)

    switch param.Sw_init
        case 1
            param = InitParams_slow(param);
            [param, psi] = cal_psi(param, psi);
            disp("after cal_psi");
            % WriteOut_data(psi, output_dir + "flux0.txt");
            [param, psi] = cal_initial_condition(param, psi);
            disp("after cal_initial_condition")
        case 2
            ReadIn_data(psi, output_dir + "flux0.txt");
            cal_initial_condition(param, psi);
        case 3
            ReadIn_data(psi, output_dir + "flux.txt");
        otherwise
            disp("Input correct number.");
    end

end

% コイルに関するデータの読み込み
function param = InitParams_slow(param)
    f = load(param.file1);
    param.coil_turn = f(1, :);
    param.coil_Ic = f(2, :) .* f(1, :);
    param.coil_z = f(3, :);
    param.coil_r = f(4, :);
    param.ncoil = size(param.coil_turn, 2);
end

function ReadIn_data(psi, filename)
    psi = readmatrix(filename);
end

function WriteOut_data(psi, filename)
    writematrix(psi, filename, 'Delimiter', 'tab');
end

function gs_graph(param, psi, p, Ip, mu_jt)
    %% グラフ
    r = param.rmin + param.delr * (0:param.Nr - 1);
    z = param.zmin + param.delz * (0:param.Nz - 1);
    [rr zz] = meshgrid(r, z);
    v = linspace(-20, 20, 21);
    fh = figure;
    set(fh, 'position', [50 50 800 600]);

    ah = subplot(1, 4, 1);
    contour(rr, zz, psi * 1000, v, 'r')

    if param.Nr == param.Nr_orig
        hold on
        psi_orig = load(comparison_C);
        contour(rr, zz, psi_orig.flux * 1000, v, 'r')
        axis equal
        legend('Matlab', 'C')
    end

    title('flux');

    v = linspace(0, 500, 21);
    ah = subplot(1, 4, 2);
    contour(rr, zz, p, v, 'k')

    if param.Nr == param.Nr_orig
        hold on
        contour(rr, zz, psi_orig.p, v, 'r')
        axis equal
        legend('Matlab', 'C')
    end

    title('presure');

    v = linspace(0.078, 0.082, 21);
    ah = subplot(1, 4, 3);
    contour(rr, zz, Ip, v, 'k')

    if param.Nr == param.Nr_orig
        hold on
        contour(rr, zz, psi_orig.rbt, v, 'r')
        axis equal
        legend('Matlab', 'C')
    end

    title('Ip');

    v = linspace(-0.5, 0, 21);
    ah = subplot(1, 4, 4);
    contour(rr, zz, mu_jt, v, 'k')

    if param.Nr == param.Nr_orig
        hold on
        contour(rr, zz, psi_orig.mujt, v, 'r')
        axis equal
        legend('Matlab', 'C')
    end

    title('mu_jt');

end

function param = define_params(param, Nz, Nr)
    param.Nz = Nz;
    param.Nr = Nr;

    % メッシュ1ステップの長さ
    param.delr = (param.rmax - param.rmin) / (Nr - 1);
    param.delz = (param.zmax - param.zmin) / (Nz - 1);

    % zに応じて格子点がいくつあるかrnumに格納
    param.rnum = zeros(1, Nz);

    for i = 1:Nz

        if ((param.zmin + (i - 1) * param.delz) < param.z1 || param.z4 < (param.zmin + (i - 1) * param.delz))
            param.rnum(1, i) = ceil((Nr - 1) * (param.rmirror - param.rmin) / (param.rmax - param.rmin));
        elseif (param.z2 < (param.zmin + (i - 1) * param.delz) && (param.zmin + (i - 1) * param.delz) < param.z3)
            param.rnum(1, i) = ceil(Nr);
        elseif (param.z1 < (param.zmin + (i - 1) * param.delz) && (param.zmin + (i - 1) * param.delz) < param.z2)
            param.rnum(1, i) = ceil((Nr - 1) * ((param.rmirror - param.rmin) / (param.rmax - param.rmin) + ...
                (param.rmax - param.rmirror) / (param.rmax - param.rmin) * (param.zmin + (i - 1) * param.delz - param.z1) / (param.z2 - param.z1)));
        elseif (param.z3 < (param.zmin + (i - 1) * param.delz) && (param.zmin + (i - 1) * param.delz) < param.z4)
            param.rnum(1, i) = ceil((Nr - 1) * (1 - (param.rmax - param.rmirror) / (param.rmax - param.rmin) * ...
                (param.zmin + (i - 1) * param.delz - param.z3) / (param.z4 - param.z3)));
        end

    end

    %% DEFINE REGIONS
    param.psi_inside = zeros(Nz, Nr);

    for i = 2:Nz - 1
        param.psi_inside(i, 2:param.rnum(1, i) - 1) = 1;
    end

    param.psi_round = zeros(Nz, Nr);

    for i = 1:Nz
        param.psi_round(i, 1:param.rnum(1, i)) = 1;
    end

    param.psi_round = param.psi_round - param.psi_inside;

    for i = 2:Nr

        for j = 1:Nz - 1

            if (param.psi_inside(j, i) == 0 && param.psi_inside(j + 1, i) == 1)
                param.psi_round(j, i) = 1;
            elseif (param.psi_inside(j, i) == 1 && param.psi_inside(j + 1, i) == 0)
                param.psi_round(j + 1, i) = 1;
            end

        end

    end

    param.psi_whole = param.psi_inside + param.psi_round;

    % リミターの位置
    if (param.rlimiter < param.rmin || param.rmax < param.rlimiter || param.zlimiter < param.zmin || param.zmax < param.zlimiter)
        % メッシュ外にリミターがあるなら(0,0)にする
        param.r_limiterpos = 0;
        param.z_limiterpos = 0;
    else
        % メッシュ内にリミターがあるなら最も近い格子点を探し，そこを指定する
        mtrxz = (0:1:Nz - 1)' .* param.delz + param.zmin;
        mtrxr = (0:1:Nr - 1) .* param.delr + param.rmin;
        dot = (mtrxz - param.zlimiter).^2 * ones(1, Nr) + ones(Nz, 1) * (mtrxr - param.rlimiter).^2;
        % 以下は2次元配列の最小値の座標を求めていて，Cに最小値を代入している
        [C, I] = min(dot(:));
        % Cは0オリジンだから値が1ずれる．
        [param.z_limiterpos, param.r_limiterpos] = ind2sub(size(dot), I);
    end

end

function [param, psi, mu_jt, p, Ip] = sub_routine_init(param, divz, divr)
    param.Nz = round(param.Nz_orig / divz);
    param.Nr = round(param.Nr_orig / divr);
    param = define_params(param, param.Nz, param.Nr);
    psi = zeros(param.Nz, param.Nr);
    mu_jt = zeros(param.Nz, param.Nr);
    p = zeros(param.Nz, param.Nr);
    Ip = zeros(param.Nz, param.Nr);
    [param, psi] = cal_initial_flux(param, psi);
end

function [param, psi, mu_jt, p, Ip] = sub_routine(param, divz, divr, psi, mu_jt, p, Ip)
    %% サブのルーチン　収束を早める
    param.Nz = round(param.Nz_orig / divz);
    param.Nr = round(param.Nr_orig / divr);
    param = define_params(param, param.Nz, param.Nr);
    disp(['    sub routine param.Nz:' num2str(param.Nz) ' param.Nr:' num2str(param.Nr)])

    psi = imresize(psi, [param.Nz, param.Nr]);
    psi = psi_round_revive(psi, param);
    mu_jt = imresize(mu_jt, [param.Nz, param.Nr]);
    p = imresize(p, [param.Nz, param.Nr]);
    Ip = imresize(Ip, [param.Nz, param.Nr]);

    for i = 1:100000
        [psi_bar, mu_jt, p, Ip] = cal_jt(param, psi, mu_jt, p, Ip);
        %     mu_jt = addPFcurrent(param,psi,mu_jt);
        format long
        [error, psi] = cal_flux_C(param, psi, mu_jt);

        if (rem(i, 100) == 0)
            disp(['    ' num2str(i) ':' num2str(error)])
        end

        if error < param.errormax
            param.routine_num = param.routine_num + i;
            disp(['  Sub routine finished at ' num2str(param.routine_num)]);
            break
        end

    end

end

% psi_roundを復活させる関数
function psi = psi_round_revive(psi, param)
    psi_b = psi;
    psi = zeros(param.Nz, param.Nr);
    [param, psi] = cal_psi(param, psi);
    psi = psi .* param.psi_round + psi_b .* (1 - param.psi_round);
end
