function run_diffusion_batch(inputs, previewFolder, analyst)

mainDir = '\\rds6.chmccorp.cchmc.org\PulMed-54\CPIR_Images_Database';
WoodsDir = '\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images';

nSubjects = height(inputs);

clc;

for i = 1:nSubjects

    fprintf('Processing subject %d of %d\n', i, nSubjects);

    currentInput = inputs(i,:);

    AgeValue      = currentInput.Age;
    SexValue      = currentInput.Sex;
    DiseaseValue  = currentInput.Disease;
    NoteValue     = currentInput.Note;
    ImageQValue   = currentInput.Quality;
    DiffFileValue = currentInput.XePath;

    if ismissing(AgeValue)
        AgeValue = "";
    end
    if ismissing(SexValue)
        SexValue = "";
    end
    if ismissing(DiseaseValue)
        DiseaseValue = "";
    end
    if ismissing(NoteValue)
        NoteValue = "";
    end
    if ismissing(ImageQValue)
        ImageQValue = "";
    end

    MainInput = struct();

    MainInput.SubjectID        = currentInput.SubjectID;
    MainInput.Age              = AgeValue;
    MainInput.Sex              = SexValue;
    MainInput.Disease          = DiseaseValue;
    MainInput.ScanDate         = currentInput.Date;
    MainInput.denoiseXe        = 'no';
    MainInput.Analyst          = analyst;
    MainInput.Polarizer        = currentInput.Polarizer;
    MainInput.N4Bias           = 'no';
    MainInput.AnalysisMethod   = 'W.Linear';
    MainInput.AgeCor           = 'no';
    MainInput.Note             = NoteValue;
    MainInput.ProcessingNotes  = currentInput.ProcessingNotes;
    MainInput.ImageQuality     = ImageQValue;
    MainInput.MaskPath         = '';%currentInput.MaskPath;
    MainInput.SliceOrientation = 'transversal';
    MainInput.DiffAcqOrder     = 'b-value interleave';
    MainInput.PreviewFolder    = previewFolder;
    MainInput.AnalysisStatus   = currentInput.AnalysisStatus;

    MainInput.SequenceType     = currentInput.Sequence;
    MainInput.SequenceName     = currentInput.SequenceName;
    sequence = string(currentInput.Sequence);    
    SequenceName = string(currentInput.SequenceName); 
    if contains(sequence, "cartesian", 'IgnoreCase', true) || contains(sequence, "RECTILINEAR", 'IgnoreCase', true)
        MainInput.SequenceType = '2D GRE';
    elseif contains(SequenceName, "vent_2d_vardens_sos", 'IgnoreCase', true)
        MainInput.SequenceType = '2D Spiral';  
    elseif contains(SequenceName, "CPIR_VENT_FLORET", 'IgnoreCase', true)
        MainInput.SequenceType = '3D FLORET';         
    elseif contains(sequence, "spiral", 'IgnoreCase', true)
        MainInput.SequenceType = '2D Spiral';
    else
        MainInput.SequenceType = char(sequence);
    end

    if contains(currentInput.Scanner, "Philips", 'IgnoreCase', true)
        MainInput.Scanner = 'Philips';     
    elseif contains(currentInput.Scanner, "GE", 'IgnoreCase', true)
        MainInput.Scanner = 'GE';
    else
        MainInput.Scanner = '';
    end    
    % Diff file
    if ismissing(DiffFileValue) || ...
            strlength(strtrim(string(DiffFileValue))) == 0
        MainInput.diff_file = '';
    else
        DiffFileValue = char(string(DiffFileValue));
        if contains(DiffFileValue, mainDir, 'IgnoreCase', true) || ...
                contains(DiffFileValue, WoodsDir, 'IgnoreCase', true)
            MainInput.diff_file = DiffFileValue;
        else
            MainInput.diff_file = fullfile(mainDir, DiffFileValue);
        end
    end

    MainInput.diff_file = char(MainInput.diff_file);

    % Get folder path from file
    [filePath,~,~] = fileparts(MainInput.diff_file);

    analysisFolder = fullfile(filePath,'Diffusion_Analysis');

    if exist(analysisFolder,'dir')
        rmdir(analysisFolder,'s');
    end

    MainInput.ScannerSoftware = char(string(currentInput.Software));

    if contains(MainInput.ScannerSoftware, "5.9")
        MainInput.ScannerSoftware = '5.9.0';
    elseif contains(MainInput.ScannerSoftware, "5.6")
        MainInput.ScannerSoftware = '5.6.1';
    elseif contains(MainInput.ScannerSoftware, "5.3")
        MainInput.ScannerSoftware = '5.3.1';
    elseif contains(MainInput.ScannerSoftware, "5.1")
        MainInput.ScannerSoftware = '5.1.7';
    elseif contains(MainInput.ScannerSoftware, "3.2")
        MainInput.ScannerSoftware = '3.2.3';
    end

    % ============================================================
    % Reconstruction type
    % ============================================================
    ACQ_Type = char(string(currentInput.SequenceName));

    if contains(ACQ_Type, "cpir_diffusion_4bs", 'IgnoreCase', true)
        num_extra_attr1 = 4;
    elseif contains(ACQ_Type, "cpir_diffusion_3bs", 'IgnoreCase', true)
        num_extra_attr1 = 3;
    elseif contains(ACQ_Type, "cpir_diffusion", 'IgnoreCase', true)
        num_extra_attr1 = 3;
    elseif contains(ACQ_Type, "hng_xe_diffusion_29jul2016", 'IgnoreCase', true)
        num_extra_attr1 = 5;
    else
        num_extra_attr1 = NaN;
    end
    % num_extra_attr1 = 4; % hard code until i figure out how to get this with the current webapp
    MainInput.sernum = num2str(currentInput.AnalysisXeSeries);
    MainInput.XeSeries = num2str(currentInput.XeSeries);

    reconValue = string(currentInput.Recon);

    if strcmpi(reconValue, "Online")
        MainInput.ReconType = 'online';
    elseif strcmpi(reconValue, "Offline")
        MainInput.ReconType = 'offline';
    end

    analysisversion = 'diff_v100';
    MainInput.analysisversion = analysisversion;

    % Format subject number with leading zeros
    subnum = currentInput.Subject;
    MainInput.subnum = subnum;

    if isnumeric(subnum)
        subnum = num2str(subnum, '%04d');
    else
        subnum = char(string(subnum));

        subnumNumeric = str2double(subnum);

        if ~isnan(subnumNumeric)
            subnum = num2str(subnumNumeric, '%04d');
        end
    end

    MainInput.num_b_values = num_extra_attr1;
    MainInput.Nbvalues     = num_extra_attr1;

    MainInput.EncodingOrder = 'linear'; % confirm this

    % ============================================================
    % Session number from scan date
    % ============================================================
    scanDate = currentInput.Date;

    if ~isdatetime(scanDate)
        scanDate = datetime(string(scanDate));
    end

    sesnum = datestr(scanDate, 'yyyymmdd');
    MainInput.ScanDate = sesnum;
    
    % Change to analysis folder and then run
    analysisfolder = fullfile( ...
        mainDir, ...
        char(string(currentInput.Study)), ...
        'analysis', ...
        analysisversion, ...
        ['sub-', subnum], ...
        ['ses-', sesnum], ...
        ['ser-', char(string(MainInput.sernum))]);

    MainInput.analysisfolder = analysisfolder;
    
    diff_analysis_folder = ...
        fullfile(analysisfolder,'Diffusion_Analysis');

    MainInput.diff_analysis_folder = ...
        diff_analysis_folder;

    if ~exist(diff_analysis_folder, 'dir')
        mkdir(diff_analysis_folder);
    end

    MainInput.OutputPath = analysisfolder;

    cd(MainInput.analysisfolder);

    % ============================================================
    % Run pipeline
    % ============================================================

    runMode = currentInput.RunMode;

    if iscell(runMode)
        runMode = runMode{1};
    end

    if isstring(runMode) || ischar(runMode)
        runMode = str2double(string(runMode));
    end

    if isempty(runMode) || ismissing(runMode) || isnan(runMode)
        runMode = 1;
    end

    if isempty(MainInput.diff_file)
        DiffusionFunctions.CCHMC_Db_Diff_Pipeline_NoData(MainInput)
    else
        if runMode == 1 && ~strcmp(string(ImageQValue), '1-Failed')
            DiffusionFunctions.CCHMC_Db_Diff_Pipeline(MainInput);
        elseif strcmp(string(ImageQValue), '1-Failed')
            DiffusionFunctions.CCHMC_Db_Diff_Pipeline_NoData(MainInput)
        elseif runMode == 2
            DiffusionFunctions.CCHMC_Db_Diff_Pipeline_rerun(MainInput);
        end

    end

end

end