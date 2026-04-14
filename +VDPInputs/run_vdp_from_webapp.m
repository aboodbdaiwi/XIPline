function run_vdp_from_webapp(excelFile, includeBackupRows)
%RUN_VDP_FROM_WEBAPP Run ventilation pipeline from web-app Excel export.
%
% Usage:
%   VDPInputs.run_vdp_from_webapp
%   VDPInputs.run_vdp_from_webapp(excelFile)
%
% Behavior:
%   - If excelFile is not provided, search Downloads for:
%       VDP_Inputs.xlsx
%       VDP_Inputs (1).xlsx
%       VDP_Inputs (2).xlsx
%       ...
%   - If exactly one matching file is found, use it.
%   - If zero or multiple matches are found, fall back to uigetfile.
%
% Notes:
%   - This entrypoint is intended to isolate web-app input handling from
%     the legacy VentilationFunctions main script.
%   - It preserves the existing downstream calls:
%       CCHMC_Db_Vent_Pipeline
%       CCHMC_Db_Vent_Pipeline_rerun
%       CCHMC_Db_Vent_Pipeline_rerun_registration
%       CCHMC_Db_Vent_Pipeline_NoData

    clc;

    % -------- Config --------
    cfg = VDPInputs.getPipelineConfig();

    if nargin < 2 || isempty(includeBackupRows)
        includeBackupRows = false;   % default: include backup rows
    end


    % -------- Resolve Excel file --------
    if nargin < 1 || isempty(excelFile) || strlength(string(excelFile)) == 0
        excelFile = iResolveExcelFile();
    else
        excelFile = char(string(excelFile));
    end

    fprintf('Using Excel file:\n%s\n\n', excelFile);

    % -------- Load / standardize / validate --------
    T = VDPInputs.LoadTable(excelFile);
    T = VDPInputs.enforceCanonicalColumns(T);
    VDPInputs.validateStandardTable(T);

    % -------- Filter runnable rows --------
    T = filter_runnable_rows(T, cfg, includeBackupRows);

    % -------- Main loop --------
    for i = 1:height(T)
        fprintf('Processing subject %d of %d\n', i, height(T));

        try
            row = T(i,:);
            MainInput = iBuildMainInput(row, cfg);

            % Run pipeline
            if iIsMissingOrEmpty(row.Vent_path{1})
                VentilationFunctions.CCHMC_Db_Vent_Pipeline_NoData(MainInput);
            else
                switch row.RunMode{1}
                    case 0
                        VentilationFunctions.CCHMC_Db_Vent_Pipeline(MainInput);
                    case 1
                        VentilationFunctions.CCHMC_Db_Vent_Pipeline_rerun(MainInput);
                    case 3
                        VentilationFunctions.CCHMC_Db_Vent_Pipeline_rerun_registration(MainInput);
                    otherwise
                        continue
                end
            end

        catch ME
            warning('VDPInputs:RowFailed', ...
                'Row %d failed: %s', i, ME.message);
            fprintf(2, '%s\n', getReport(ME, 'basic', 'hyperlinks', 'off'));
        end
    end
end


function excelFile = iResolveExcelFile()
    % Look in Downloads for VDP_Inputs*.xlsx.
    downloadsDir = fullfile(getenv('USERPROFILE'), 'Downloads');

    if ~isfolder(downloadsDir)
        excelFile = iPromptForExcelFile();
        return
    end

    listing = dir(fullfile(downloadsDir, 'VDP_Inputs*.xlsx'));

    % Keep only:
    %   VDP_Inputs.xlsx
    %   VDP_Inputs (N).xlsx
    keep = false(numel(listing), 1);
    for k = 1:numel(listing)
        keep(k) = ~isempty(regexp(listing(k).name, '^VDP_Inputs( \(\d+\))?\.xlsx$', 'once'));
    end
    listing = listing(keep);

    if isempty(listing)
        excelFile = iPromptForExcelFile();
        return
    end

    % Choose the most recently modified file
    [~, idx] = max([listing.datenum]);
    latest = listing(idx);
    listing(idx) = [];

    excelFile = fullfile(latest.folder, latest.name);
    
    fprintf('Multiple matching files found. Using most recent:\n%s\n\n', excelFile);

    fprintf('Deleting old VDP input files:\n');
    for k = 1:numel(listing)
        fpath = fullfile(listing(k).folder, listing(k).name);
        try
            delete(fpath);
        catch ME
            warning('VDPInputs:DeleteFailed', ...
                'Could not delete file: %s\nReason: %s', ...
                fpath, ME.message);
        end
    end
end


function excelFile = iPromptForExcelFile()
    [f, p] = uigetfile({'*.xlsx', 'Excel Files (*.xlsx)'}, ...
        'Select VDP input Excel file');

    if isequal(f, 0)
        error('VDPInputs:NoExcelSelected', ...
            'No Excel file was selected.');
    end

    excelFile = fullfile(p, f);
end


function MainInput = iBuildMainInput(row, cfg)
    MainInput = struct();

    % Core metadata
    MainInput.SubjectID        = iCellDefault(row.SubjectID{1}, "");
    MainInput.Age              = iCellDefault(row.AGE{1}, "");
    MainInput.Sex              = iCellDefault(row.Sex{1}, "");
    MainInput.Disease          = iCellDefault(row.DiseaseType{1}, "");
    MainInput.ScanDate         = iCellDefault(row.ScanDate{1}, "");
    MainInput.Scanner          = 'Philips-3T';
    MainInput.ScannerSoftware  = iCellDefault(row.ScannerSoftware{1}, "");
    MainInput.SequenceType     = iCellDefault(row.TRAVERSAL_GEO{1}, "");
    MainInput.denoiseXe        = 'no';
    MainInput.Analyst          = cfg.analyzer;

    % Voxel / orientation
    vx = iCellDefault(row.PIXEL_SPACING_X{1}, "");
    vy = iCellDefault(row.PIXEL_SPACING_Y{1}, "");
    vz = iCellDefault(row.SLICE_THICKNESS{1}, "");

    MainInput.VoxelSize        = sprintf('[%s,%s,%s]', ...
                                  num2str(vx), num2str(vy), num2str(vz));
    MainInput.PIXEL_SPACING_X  = vx;
    MainInput.PIXEL_SPACING_Y  = vy;
    MainInput.SLICE_THICKNESS  = vz;
    MainInput.SliceOrientation = iCellDefault(row.SliceOrientation{1}, "");

    % Optional fields
    protonSer   = iCellDefault(row.ProtonSeriesNumber{1}, "");
    anatPathRaw = iCellDefault(row.PROT_FILEPATH{1}, "");
    imageQ      = iCellDefault(row.ImageQuality{1}, "");
    noteVal     = iCellDefault(row.Notes{1}, "");

    % Vent / anat file paths
    ventPathRaw = iCellDefault(row.Vent_path{1}, "");
    if iIsMissingOrEmpty(ventPathRaw)
        MainInput.vent_file = '';
        MainInput.anat_file = '';
        MainInput.ImageQuality = '1-Failed';
    else
        MainInput.vent_file = iResolveDataPath(ventPathRaw, cfg);
        MainInput.anat_file = iResolveDataPath(anatPathRaw, cfg);
        MainInput.ImageQuality = imageQ;
    end

    MainInput.vent_file     = char(string(MainInput.vent_file));
    MainInput.anat_file     = char(string(MainInput.anat_file));

    % Recon type
    [~, ~, ext] = fileparts(MainInput.vent_file);
    if strcmpi(ext, '.dcm')
        MainInput.ReconType = 'online';
    elseif strcmpi(ext, '.data')
        MainInput.ReconType = 'offline';
    else
        MainInput.ReconType = '';
    end

    MainInput.analysisversion = cfg.analysisVersion;

    % Subject number with leading zeros
    subnum = row.Subject{1};
    if isnumeric(subnum)
        subnum = num2str(subnum, '%04d');
    else
        subnum = char(string(subnum));
    end

    MainInput.SCAN_NUM         = row.Scan_num;
    MainInput.xe_sernum        = row.XeSeriesNumber{1};
    MainInput.proton_sernum    = protonSer;
    MainInput.Note             = noteVal;
    MainInput.freqoffset       = iCellDefault(row.FrequencyOffsetHz{1}, "");

    % Analysis folder
    studyVal = iCellDefault(row.Study{1}, "");
    sesNum   = row.ScanDate{1};

    MainInput.analysisFolder = fullfile( ...
        cfg.mainDir, ...
        studyVal, ...
        'analysis', ...
        cfg.analysisVersion, ...
        ['sub-', subnum], ...
        ['ses-', num2str(sesNum)], ...
        ['ser-', num2str(MainInput.xe_sernum)]);

    if ~exist(MainInput.analysisFolder, 'dir')
        mkdir(MainInput.analysisFolder);
    end

    cd(MainInput.analysisFolder);
end


function p = iResolveDataPath(p, cfg)
    % Normalize mixed path formats from new web-app table.

    p = char(string(p));
    p = strtrim(p);

    if isempty(p)
        return
    end

    % Normalize separators
    p = strrep(p, '/', '\');

    % Fix leading "/\\" case
    if startsWith(p, '/\\')
        p = p(2:end);
    end

    % Already absolute UNC
    if startsWith(p, '\\')
        return
    end

    % Already drive-letter path
    if ~isempty(regexp(p, '^[A-Za-z]:\\', 'once'))
        return
    end

    % Looks like a network path missing UNC prefix
    if contains(p, 'PulMed-54', 'IgnoreCase', true) || ...
       contains(p, 'PulMed-35', 'IgnoreCase', true) || ...
       startsWith(p, 'rds', 'IgnoreCase', true)
        p = ['\\' p];
        return
    end

    % If it already looks like mainDir tail or WoodsDir, preserve intent
    if contains(fullfile(p), cfg.mainDir(17:end), 'IgnoreCase', true) || ...
       contains(p, cfg.WoodsDir, 'IgnoreCase', true)
        return
    end

    % Otherwise treat as relative to mainDir
    p = fullfile(cfg.mainDir, p);
end


function tf = iIsMissingOrEmpty(x)
    if ismissing(x)
        tf = true;
        return
    end

    if isstring(x) || ischar(x)
        tf = strlength(string(x)) == 0;
        return
    end

    tf = isempty(x);
end


function v = iCellDefault(v, defaultValue)
    if ismissing(v)
        v = defaultValue;
        return
    end

    if (isstring(v) || ischar(v)) && strlength(string(v)) == 0
        v = defaultValue;
    end
end