function CONFIG = loadConfigure(ConfigPath)
    A = load(ConfigPath)
    CONFIG.DataType = A(1, 1);
    CONFIG.DevideFlux = A(2, 1);
    CONFIG.DetermineCCSZPos = A(3, 1);
    CONFIG.DetermineCCSRPos = A(4, 1);
    CONFIG.ShowFig = A(5, 1);

end
