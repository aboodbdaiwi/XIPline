function T_out = filter_runnable_rows(T, cfg, includeBackupRows)

    runMode = cell2mat(T.RunMode);
    T = T(~isnan(runMode) & runMode ~= 2, :);

    polarizer = string(T.Polarizer);

    isMyPrimary    = polarizer == string(cfg.analyzerInfo.primaryPolarizer);
    isOtherPrimary = polarizer == string(cfg.analyzerInfo.otherPolarizer);
    isBackup       = ~(isMyPrimary | isOtherPrimary);
    
    if includeBackupRows
        T_out = T(isMyPrimary | isBackup, :);
    else
        T_out = T(isMyPrimary, :);
    end
    
    fprintf('Runnable rows: %d\n\n', height(T_out));

end