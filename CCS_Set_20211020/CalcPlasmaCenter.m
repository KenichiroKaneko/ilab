function CCS_R = CalcPlasmaCenter(PARAM, CONFIG, CCS_Z)
    % インボード側の磁束の極大値から求めたCCS_Zと
    % GS方程式から求めたプラズマ中心の時間発展の情報
    % を利用して、CCS_Rを決定する
    % PlasmaCenterRZ.matを読み込む必要がある
    % CCScenterRZは４列のデータ[R Z r z]から成る
    % RZはメッシュ点、rzはメートル単位のUTST内の座標
    % 第1象限のみから成る配列なので最初に拡張している。

    CCScenterRZ = load("CCS_temporary/PlasmaCenterRZ").CCScenterRZ;
    r = CCScenterRZ(:, 3);
    z = CCScenterRZ(:, 4);
    num = length(z);
    z = [flipud(z); -z]; r = [flipud(r); r];
    z(num) = []; r(num) = [];

    for i = 1:PARAM.CCS
        CCS_R(i) = spline(z, r, CCS_Z(i));
    end

    if CONFIG.ShowFig
        figure()
        title("プラズマ中心の時間発展とスプライン補完によるCCS面中心の決定")
        f = fit(z, r, "smoothingspline");
        fnplt(f.p, "-k", 1)
        hold on

        for i = 1:PARAM.CCS
            CCS_R(i) = spline(z, r, CCS_Z(i));
            plot(z, r, "om", CCS_Z(i), CCS_R(i), "*r")
        end

        view(90, 90);
        hold off
    end

end
