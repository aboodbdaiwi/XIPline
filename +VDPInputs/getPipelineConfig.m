function cfg = getPipelineConfig()

    cfg = struct();

    cfg.mainDir         = '\\Rds6.cchmc.org\PulMed-54\CPIR_Images_Database';
    cfg.WoodsDir        = '\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images';
    cfg.analysisVersion = 'vent_v100';
    
    cfg.analyzerInfo    = VDPInputs.detectAnalyzer();
    cfg.analyzer        = cfg.analyzerInfo.code;

end