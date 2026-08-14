function run_ventilation_batch(inputs, previewFolder, analyst)

mainDir = '\\rds6.chmccorp.cchmc.org\PulMed-54\CPIR_Images_Database';
WoodsDir = '\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images';

% ================================================================
% Inputs now come directly from the app
% inputs = app.BatchProcessing.filtered_inputs
% ================================================================

nSubjects = height(inputs);

clc;

for i = 1:nSubjects

    fprintf('Processing subject %d of %d\n', i, nSubjects);

    % ============================================================
    % Get current row
    % ============================================================
    currentInput = inputs(i,:);

    % ============================================================
    % Handle missing values
    % ============================================================
    AgeValue = currentInput.Age;
    SexValue = currentInput.Sex;
    HSerNumValue = currentInput.HSeries;
    AnatFileValue = currentInput.HPath;
    ImageQValue = currentInput.Quality;
    NoteValue = currentInput.Note;
    

    if ismissing(AgeValue)
        AgeValue = "";
    end

    if ismissing(SexValue)
        SexValue = "";
    end

    if ismissing(HSerNumValue)
        HSerNumValue = "";
    end

    if ismissing(AnatFileValue)
        AnatFileValue = "";
    end

    if ismissing(ImageQValue)
        ImageQValue = "";
    end

    if ismissing(NoteValue)
        NoteValue = "";
    end

    % ============================================================
    % Main input structure
    % ============================================================
    MainInput = struct();

    MainInput.SubjectID        = currentInput.SubjectID;
    MainInput.Age              = AgeValue;
    MainInput.Sex              = SexValue;
    MainInput.Disease          = currentInput.Disease;
    MainInput.ScanDate         = currentInput.Date;
    MainInput.AnalysisStatus   = currentInput.AnalysisStatus;

    % Preserve previous default if Scanner is missing
    if ismissing(currentInput.Scanner) || ...
            strlength(strtrim(string(currentInput.Scanner))) == 0
        MainInput.Scanner = 'Philips-3T';
    else
        MainInput.Scanner = currentInput.Scanner;
    end

    MainInput.ScannerSoftware  = currentInput.Software;
    MainInput.SequenceType     = currentInput.Sequence;
    sequence = string(currentInput.Sequence);    
    if contains(sequence, "cartesian", 'IgnoreCase', true)
        MainInput.SequenceType = '2D GRE';
    elseif contains(sequence, "spiral", 'IgnoreCase', true)
        MainInput.SequenceType = '2D Spiral';
    else
        MainInput.SequenceType = char(sequence);
    end

    MainInput.denoiseXe        = 'no';

    % Analyst now comes from Polarizer
    MainInput.Polarizer        = currentInput.Polarizer;
    MainInput.Analyst          = analyst;

    MainInput.VoxelSize = sprintf( ...
        '[%s,%s,%s]', ...
        num2str(currentInput.PixelX), ...
        num2str(currentInput.PixelY), ...
        num2str(currentInput.SliceThickness));

    MainInput.PIXEL_SPACING_X  = currentInput.PixelX;
    MainInput.PIXEL_SPACING_Y  = currentInput.PixelY;
    MainInput.SLICE_THICKNESS  = currentInput.SliceThickness;
    MainInput.SliceOrientation = currentInput.Orientation;

    % ============================================================
    % Vent file
    % ============================================================
    VentFileValue = currentInput.XePath;

    if ismissing(VentFileValue) || ...
            strlength(strtrim(string(VentFileValue))) == 0

        MainInput.vent_file = '';
        MainInput.anat_file = '';
        ImageQValue = '1-Failed';

    else

        VentFileValue = char(string(VentFileValue));

        if contains(VentFileValue, mainDir, 'IgnoreCase', true) || ...
                contains(VentFileValue, WoodsDir, 'IgnoreCase', true)

            MainInput.vent_file = VentFileValue;

        else

            MainInput.vent_file = fullfile(mainDir, VentFileValue);

        end

        % ========================================================
        % Anat file
        % ========================================================
        if ismissing(AnatFileValue) || ...
                strlength(strtrim(string(AnatFileValue))) == 0

            MainInput.anat_file = '';

        else

            AnatFileValue = char(string(AnatFileValue));

            if contains(AnatFileValue, mainDir, 'IgnoreCase', true) || ...
                    contains(AnatFileValue, WoodsDir, 'IgnoreCase', true)

                MainInput.anat_file = AnatFileValue;

            else

                MainInput.anat_file = fullfile(mainDir, AnatFileValue);

            end

        end

    end

    MainInput.SCAN_NUM  = currentInput.Scan;
    MainInput.vent_file = char(MainInput.vent_file);
    MainInput.anat_file = char(MainInput.anat_file);

    % ============================================================
    % Recon type
    %
    % Use the new Recon column instead of file extension because
    % offline paths may end in .*
    % ============================================================
    reconValue = string(currentInput.Recon);

    if strcmpi(reconValue, "Online")
        MainInput.ReconType = 'online';

    elseif strcmpi(reconValue, "Offline")
        MainInput.ReconType = 'offline';

    else
        MainInput.ReconType = '';
    end

    analysisversion = 'vent_v100';
    MainInput.analysisversion = analysisversion;

    % ============================================================
    % Format subject number with leading zeros
    % ============================================================
    subnum = currentInput.Subject;

    if isnumeric(subnum)
        subnum = num2str(subnum, '%04d');
    else
        subnum = char(string(subnum));

        subnumNumeric = str2double(subnum);

        if ~isnan(subnumNumeric)
            subnum = num2str(subnumNumeric, '%04d');
        end
    end

    % ============================================================
    % Session number from ScanDate
    % ============================================================
    scanDate = currentInput.Date;

    if ~isdatetime(scanDate)
        scanDate = datetime(string(scanDate));
    end

    sesnum = datestr(scanDate, 'yyyymmdd');
    MainInput.ScanDate = sesnum;
    % ============================================================
    % Remaining inputs
    % ============================================================
    MainInput.xe_sernum     = currentInput.XeSeries;
    MainInput.proton_sernum = HSerNumValue;
    MainInput.ImageQuality  = ImageQValue;
    MainInput.Note          = NoteValue;
    MainInput.ProcessingNotes = currentInput.ProcessingNotes;
    MainInput.freqoffset    = currentInput.FreqOffset;

    % ============================================================
    % Preview folder passed from app
    % ============================================================
    MainInput.PreviewFolder = previewFolder;

    % ============================================================
    % Analysis folder
    % ============================================================
    MainInput.analysisFolder = fullfile( ...
        mainDir, ...
        char(string(currentInput.Study)), ...
        'analysis', ...
        analysisversion, ...
        ['sub-', subnum], ...
        ['ses-', sesnum], ...
        ['ser-', num2str(MainInput.xe_sernum)]);

    if ~exist(MainInput.analysisFolder, 'dir')
        mkdir(MainInput.analysisFolder);
    end

    cd(MainInput.analysisFolder);

    % ============================================================
    % Run pipeline
    % ============================================================
    if isempty(MainInput.vent_file)

        VentilationFunctions.CCHMC_Db_Vent_Pipeline_NoData(MainInput)

    else

        runMode = currentInput.RunMode;

        % Convert RunMode safely if it came in as string/text
        if iscell(runMode)
            runMode = runMode{1};
        end

        if isstring(runMode) || ischar(runMode)
            runMode = str2double(string(runMode));
        end

        % Missing RunMode -> default to normal run
        if isempty(runMode) || ismissing(runMode) || isnan(runMode)
            runMode = 1;
        end

        if runMode == 1
            VentilationFunctions.CCHMC_Db_Vent_Pipeline(MainInput);
        elseif runMode == 2
            VentilationFunctions.CCHMC_Db_Vent_Pipeline_rerun(MainInput);
        elseif runMode == 3
            VentilationFunctions.CCHMC_Db_Vent_Pipeline_rerun_registration(MainInput);
        end

    end

end

end