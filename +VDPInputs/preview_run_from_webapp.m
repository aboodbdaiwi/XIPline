function preview_run_from_webapp(excelFile, includeBackupRows)

    if nargin < 1 || isempty(excelFile)
        excelFile = VDPInputs.resolve_excel_file();
    end
    if nargin < 2 || isempty(includeBackupRows)
        includeBackupRows = false;
    end

    cfg = VDPInputs.getPipelineConfig();

    T = VDPInputs.LoadTable(excelFile);
    T = VDPInputs.enforceCanonicalColumns(T);
    VDPInputs.validateStandardTable(T);

    Tpreview = VDPInputs.filter_runnable_rows(T, cfg, includeBackupRows);

    fprintf('Excel file:\n%s\n\n', excelFile);
    fprintf('Detected analyzer: %s\n', cfg.analyzerInfo.name);
    fprintf('Primary polarizer: %s\n', cfg.analyzerInfo.primaryPolarizer);
    fprintf('Include backup rows: %d\n', includeBackupRows);
    fprintf('Rows to process: %d\n\n', height(Tpreview));

    disp(unique(string(Tpreview.Polarizer)));

    assignin('base', 'VDPPreviewTable', Tpreview);
    openvar('VDPPreviewTable');

    

end