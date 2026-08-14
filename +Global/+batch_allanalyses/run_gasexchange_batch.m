function run_gasexchange_batch(inputs, previewFolder, analyst)

mainDir = '\\rds6.chmccorp.cchmc.org\PulMed-54\CPIR_Images_Database';
WoodsDir = '\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images';

nSubjects = height(inputs);

clc;

for i = 1:nSubjects

    fprintf('Processing subject %d of %d\n', i, nSubjects);

    currentInput = inputs(i,:);

    % ============================================================
    % Handle missing values
    % ============================================================
    AgeValue     = currentInput.Age;
    SexValue     = currentInput.Sex;
    DiseaseValue = currentInput.Disease;
    NoteValue    = currentInput.Note;
    % CalFileValue = currentInput.CalPath;
    GxFileValue  = currentInput.XePath;
    AnatFileValue = currentInput.HPath;

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

    % if ismissing(CalFileValue)
    %     disp("no calibration .data file detected")
    %     CalFileValue = "";
    % end

    % ============================================================
    % Main input structure
    % ============================================================
    MainInput = struct();

    MainInput.SubjectID       = currentInput.SubjectID;
    MainInput.Age             = AgeValue;
    MainInput.Sex             = SexValue;
    MainInput.Disease         = DiseaseValue;
    MainInput.ScanDate        = currentInput.Date;
    MainInput.Scanner         = currentInput.Scanner;
    MainInput.AnalysisVersion = 'CCHMC';   %currentInput.Software;
    MainInput.ScannerSoftware = char(string(currentInput.Software));
    if contains(MainInput.ScannerSoftware,'-')
        MainInput.ScannerSoftware = strrep(MainInput.ScannerSoftware,'-','.');
    end    
    MainInput.SequenceType    = '3D-Radial';
    MainInput.denoiseXe       = 'no';
    MainInput.Analyst         = analyst;
    MainInput.Polarizer       = currentInput.Polarizer;
    MainInput.N4Bias          = 'yes';
    MainInput.AnalysisMethod  = '1-Point Dixon';
    MainInput.AgeCor          = 'no';
    MainInput.PreviewFolder   = previewFolder;
    MainInput.AnalysisStatus   = currentInput.AnalysisStatus;

    % ============================================================
    % Gas-exchange file
    % ============================================================
    if ismissing(GxFileValue) || ...
            strlength(strtrim(string(GxFileValue))) == 0

        MainInput.gx_file      = '';
        MainInput.anat_file    = '';
        MainInput.cal_file     = '';
        MainInput.ImageQuality = '1-Failed';

        scanDate = currentInput.Date;

        if ~isdatetime(scanDate)
            scanDate = datetime(string(scanDate));
        end

        scandate = datestr(scanDate, 'yyyymmdd');

        timeStr = '';

        MainInput.scandate = scandate;
        MainInput.timeStr  = timeStr;

        NoDataflag = 1;

    else

        MainInput.gx_file   = char(string(GxFileValue));
        MainInput.anat_file = char(string(AnatFileValue));

        % Find latest calibration .sin file and convert to corresponding .data file
        gx_file = fullfile(mainDir,MainInput.gx_file);
        gxPath = char(gx_file);
        gxFolder = fileparts(gxPath);
        sessionFolder = fileparts(gxFolder);
        calFolder = fullfile(sessionFolder,'cal');
        if ~isfolder(calFolder)
            error('Calibration folder not found: %s',calFolder);
        end
        calSinFiles = dir(fullfile(calFolder,'*.sin'));
        if isempty(calSinFiles)
            error('No calibration .sin files found in: %s',calFolder);
        end
        calTimes = NaT(numel(calSinFiles),1);
        for i = 1:numel(calSinFiles)
            sinFile = fullfile(calSinFiles(i).folder,calSinFiles(i).name);
            txt = fileread(sinFile);
            tok = regexp(txt,'start_scan_date_time\s*:\s*(\d{1,2}-[A-Za-z]{3}-\d{4})\s+(\d{2}:\d{2}:\d{2}(?:\.\d+)?)','tokens','once');
            if ~isempty(tok)
                try
                    calTimes(i) = datetime([tok{1} ' ' tok{2}],'InputFormat','dd-MMM-yyyy HH:mm:ss.SS','Locale','en_US');
                catch
                    calTimes(i) = datetime([tok{1} ' ' tok{2}],'InputFormat','dd-MMM-yyyy HH:mm:ss','Locale','en_US');
                end
            end
        end
        validFiles = ~isnat(calTimes);
        if ~any(validFiles)
            error('Could not find start_scan_date_time in any calibration .sin file.');
        end
        validIdx = find(validFiles);
        [~,idx] = max(calTimes(validFiles));
        latestIdx = validIdx(idx);
        latestCalSin = fullfile(calSinFiles(latestIdx).folder,calSinFiles(latestIdx).name);
        [calPath,calName,~] = fileparts(latestCalSin);
        MainInput.cal_file = fullfile(calPath,[calName '.data']);
        fprintf('Latest calibration time: %s\n',datestr(calTimes(latestIdx),'dd-mmm-yyyy HH:MM:SS.FFF'));
        fprintf('Calibration file selected:\n%s\n',MainInput.cal_file);
        % % Remove main directory from calibration path
        % calFile = char(MainInput.cal_file);
        % mainPath = char(mainDir);
        % if startsWith(calFile,mainPath,'IgnoreCase',true)
        %     calFile = calFile(length(mainPath)+1:end);
        % end
        % if ~startsWith(calFile,filesep)
        %     calFile = [filesep calFile];
        % end
        % MainInput.cal_file = calFile;

        if contains(MainInput.gx_file, mainDir, 'IgnoreCase', true) || ...
                contains(MainInput.gx_file, WoodsDir, 'IgnoreCase', true)
            MainInput.gx_file = char(string(GxFileValue));
        else
            MainInput.gx_file = fullfile( ...
                mainDir, char(string(GxFileValue)));
        end

        % ========================================================
        % Anat file
        % ========================================================
        if ~ismissing(AnatFileValue) && ...
                strlength(strtrim(string(AnatFileValue))) > 0

            if contains(char(string(AnatFileValue)), mainDir, 'IgnoreCase', true) || ...
                    contains(char(string(AnatFileValue)), WoodsDir, 'IgnoreCase', true)

                MainInput.anat_file = char(string(AnatFileValue));
            else

                MainInput.anat_file = fullfile( ...
                    mainDir, char(string(AnatFileValue)));
            end
        else
            MainInput.anat_file = '';
        end

        % Find scan date/time from .sin file
        [folder,~,~] = fileparts(MainInput.gx_file);
        sinfiles = dir(fullfile(folder,'**','*.sin*'));
        sinfiles = sinfiles(~[sinfiles.isdir]);
        names = {sinfiles.name};
        keep = ~contains(names,'CoilSurveyScan','IgnoreCase',true) & ~contains(names,'Mask','IgnoreCase',true);
        sinfiles = sinfiles(keep);
        if isempty(sinfiles)
            error('No eligible .sin files found after filtering.');
        end
        scanTimes = NaT(numel(sinfiles),1);
        for i = 1:numel(sinfiles)
            sinPath = fullfile(sinfiles(i).folder,sinfiles(i).name);
            txt = fileread(sinPath);
            tok = regexp(txt,'start_scan_date_time\s*:\s*(\d{1,2}-[A-Za-z]{3}-\d{4})\s+(\d{2}:\d{2}:\d{2}(?:\.\d+)?)','tokens','once');
            if ~isempty(tok)
                try
                    scanTimes(i) = datetime([tok{1} ' ' tok{2}],'InputFormat','dd-MMM-yyyy HH:mm:ss.SS','Locale','en_US');
                catch
                    try
                        scanTimes(i) = datetime([tok{1} ' ' tok{2}],'InputFormat','dd-MMM-yyyy HH:mm:ss','Locale','en_US');
                    catch
                    end
                end
            end
        end
        validFiles = ~isnat(scanTimes);
        if ~any(validFiles)
            error('Could not find start_scan_date_time in any eligible .sin file.');
        end
        validIdx = find(validFiles);
        [~,idx] = max(scanTimes(validFiles));
        latestIdx = validIdx(idx);
        scanTime = scanTimes(latestIdx);

        scandate = datestr(scanTime,'yyyymmdd');
        timeStr = datestr(scanTime,'HHMMSS');

        MainInput.scandate = scandate;
        MainInput.timeStr = timeStr;

        fprintf('Scan date/time: %s\n',datestr(scanTime,'dd-mmm-yyyy HH:MM:SS.FFF'));
        fprintf('Source SIN file: %s\n',fullfile(sinfiles(latestIdx).folder,sinfiles(latestIdx).name));

        MainInput.cal_file  = char(MainInput.cal_file);
        MainInput.gx_file   = char(MainInput.gx_file);
        MainInput.anat_file = char(MainInput.anat_file);

        NoDataflag = 0;

    end

    % ============================================================
    % Reconstruction / version
    % ============================================================
    MainInput.ReconType = 'offline';

    analysisversion = 'gx_v100';

    MainInput.analysisversion = analysisversion;
    MainInput.Note = NoteValue;
    MainInput.ProcessingNotes = currentInput.ProcessingNotes;

    % ============================================================
    % Subject number
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
    % Session number
    % ============================================================
    scanDate = currentInput.Date;

    if ~isdatetime(scanDate)
        scanDate = datetime(string(scanDate));
    end

    sesnum = datestr(scanDate, 'yyyymmdd');

    % ============================================================
    % Analysis folder
    % ============================================================
    analysisfolder = fullfile( ...
        mainDir, ...
        char(string(currentInput.Study)), ...
        'analysis', ...
        analysisversion, ...
        ['sub-', subnum], ...
        ['ses-', sesnum], ...
        ['ser-', timeStr]);

    MainInput.analysisfolder = analysisfolder;

    gx_analysis_folder = ...
        fullfile(analysisfolder, 'GasExchange_Analysis');

    MainInput.gx_analysis_folder = gx_analysis_folder;

    if ~exist(gx_analysis_folder, 'dir')
        mkdir(gx_analysis_folder);
    end

    cd(MainInput.analysisfolder);

    % ============================================================
    % For JSON
    % ============================================================
    MainInput.analysispath = fullfile( ...
        '\', ...
        char(string(currentInput.Study)), ...
        'analysis', ...
        analysisversion, ...
        ['sub-', subnum], ...
        ['ses-', sesnum], ...
        ['ser-', timeStr]);

    analysispath = MainInput.analysispath;

    % ============================================================
    % Run mode
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

    % ============================================================
    % Run pipeline
    % ============================================================
    if NoDataflag == 1
        GasExchangeFunctions.CCHMC_Db_GX_Pipeline_NoData(MainInput)
    else
        if runMode == 1
            GasExchangeFunctions.CCHMC_Db_GX_Pipeline(MainInput);
        elseif runMode == 2
            MainInput.UpdatedNote = currentInput.Note;
            MainInput.UpdatedImageQuality = currentInput.Quality;
            MainInput.UpdatedProcessingNotes = currentInput.ProcessingNotes;
            MainInput.UpdatedAnalysisStatus = currentInput.AnalysisStatus;

            GasExchangeFunctions.CCHMC_Db_GX_Pipeline_rerun( ...
                MainInput, analysispath);
        end
    end
end

end