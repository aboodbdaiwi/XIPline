clc; clear;


excelFile = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CPIR_GX_analysis\gx_main_input.xlsx';
mainDir = '\\rds6.chmccorp.cchmc.org\PulMed-54\CPIR_Images_Database';
WoodsDir = '\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images';

% === Load Excel Table ===
T = readcell(excelFile);

% Extract relevant columns into simple cell arrays
SubjectCol   = T(:,1);
ScanDateCol  = T(:,2);
StudyCol     = T(:,3);
SubNumCol    = T(:,4);
DiseaseCol   = T(:,5);
AgeCol       = T(:,6);
SexCol       = T(:,7);
ScannerCol   = T(:,8);
ScanVerCol   = T(:,9);
CalFileCol   = T(:,10);
GxFileCol    = T(:,11);
GxAnatFileCol  = T(:,12);
NoteCol      = T(:,13);
RuneCol      = T(:,14);

nSubjects = size(SexCol,1);

%% 

clc;
for i = 4%:nSubjects % always start from 2
    fprintf('Processing subject %d of %d\n', i, nSubjects);

    if ismissing(AgeCol{i})
        AgeCol{i} = "";
    end
    if ismissing(SexCol{i})
        SexCol{i} = "";
    end  
    if ismissing(DiseaseCol{i})
        DiseaseCol{i} = "";
    end                 
    if ismissing(NoteCol{i})
        NoteCol{i} = "";
    end

    MainInput = struct();
    MainInput.SubjectID        = SubjectCol{i};
    MainInput.Age              = AgeCol{i};
    MainInput.Sex              = SexCol{i}; 
    MainInput.Disease          = DiseaseCol{i}; 
    MainInput.ScanDate         = ScanDateCol{i};
    MainInput.Scanner          = ScannerCol{i}; 
    MainInput.AnalysisVersion  = ScanVerCol{i};
    MainInput.ScannerSoftware  = '5.9.0';
    MainInput.SequenceType     = '3D-Radial';
    MainInput.denoiseXe        = 'no';
    MainInput.Analyst          = 'Database';
    MainInput.N4Bias           = 'yes';
    MainInput.AnalysisMethod   = '1-Point Dixon';
    MainInput.AgeCor           = 'no';    


    % Gx file
    if ismissing(GxFileCol{i}) 
        MainInput.gx_file           = '';
        MainInput.anat_file         = '';
        MainInput.ImageQuality      = '1-Failed';
        timeStr                     = '';   
        MainInput.scandate          = ScanDateCol{i};
        MainInput.timeStr           = timeStr;    
        NoDataflag                  = 1;
    else    
        MainInput.gx_file = GxFileCol{i};
        MainInput.anat_file = GxAnatFileCol{i};
        MainInput.cal_file = char(CalFileCol{i});        
        if contains(GxFileCol{i}, mainDir, 'IgnoreCase', true) || contains(GxFileCol{i}, WoodsDir, 'IgnoreCase', true)
            MainInput.gx_file = GxFileCol{i};
        else
            MainInput.gx_file = fullfile(mainDir, GxFileCol{i});
        end
        
        % Anat file
        if contains(GxAnatFileCol{i}, mainDir, 'IgnoreCase', true) || contains(GxAnatFileCol{i}, WoodsDir, 'IgnoreCase', true)
            MainInput.anat_file = GxAnatFileCol{i};
        else
            MainInput.anat_file = fullfile(mainDir, GxAnatFileCol{i});
        end
        [folder, name, ~] = fileparts(MainInput.gx_file);

        % Get all files recursively
        sinfiles = dir(fullfile(folder, '**', '*.sin*'));
        sinfiles = sinfiles(~[sinfiles.isdir]);  % remove folders
        
        % Exclude names containing 'CoilSurveyScan' or 'Mask' (case-insensitive)
        names = {sinfiles.name};
        keep  = ~contains(names, 'CoilSurveyScan', 'IgnoreCase', true) & ...
                ~contains(names, 'Mask',          'IgnoreCase', true);
        sinfiles = sinfiles(keep);
        
        if isempty(sinfiles)
            error('No eligible .sin files found after filtering.');
        end
        % Pick the last file by modification time
        [~, order] = sort([sinfiles.datenum], 'ascend');
        lastFile   = sinfiles(order(end));
        lastPath   = fullfile(lastFile.folder, lastFile.name);
        % Extract time (2nd numeric group: yyyymmdd_HHMMSS_...)
        tok = regexp(lastFile.name, '^(\d{8})_(\d{6})_', 'tokens', 'once');
        if isempty(tok)
            error('Filename does not match expected pattern: yyyymmdd_HHMMSS_...');
        end
        scandate = tok{1};
        timeStr = tok{2};           % e.g., '132833'
        MainInput.scandate     = scandate;
        MainInput.timeStr     = timeStr; 
        MainInput.cal_file     = char(MainInput.cal_file);
        MainInput.gx_file     = char(MainInput.gx_file);
        MainInput.anat_file     = char(MainInput.anat_file);
        NoDataflag = 0;
    end

    MainInput.ReconType = 'offline';
    analysisversion = 'gx_v100';
    MainInput.analysisversion = analysisversion;
    MainInput.Note             = NoteCol{i};

    % make the new analysis folder location
    subnum = SubNumCol{i};
    if isnumeric(subnum)
        subnum = num2str(subnum, '%04d');
    end
    
    analysisfolder = fullfile(mainDir,StudyCol{i}, 'analysis', analysisversion, ...
        ['sub-', subnum], ['ses-', num2str(MainInput.ScanDate)], ['ser-', timeStr]);
    MainInput.analysisfolder = analysisfolder;

    gx_analysis_folder = fullfile(analysisfolder,'GasExchange_Analysis');
    MainInput.gx_analysis_folder = gx_analysis_folder;
    if ~exist(gx_analysis_folder, 'dir')
        mkdir(gx_analysis_folder);
    end
    cd(MainInput.analysisfolder);

    % Run pipeline
    if NoDataflag == 1
        GasExchangeFunctions.CCHMC_Db_GX_Pipeline_NoData(MainInput)
    else
        if RuneCol{i} == 0
            GasExchangeFunctions.CCHMC_Db_GX_Pipeline(MainInput);
        elseif RuneCol{i} == 1
            GasExchangeFunctions.CCHMC_Db_GX_Pipeline_rerun(MainInput);
        end
        % T{i,15} = MainInput.analysisfolder;
    end
end
% % Define output file path
% [folder, name, ~] = fileparts(excelFile);
% outFile = fullfile(folder, [name '_withpath.xlsx']);
% 
% % Replace all `missing` values in the cell array with empty strings
% for r = 1:numel(T)
%     if ismissing(T{r})
%         T{r} = '';   % or 'NA' if you prefer
%     end
% end
% 
% % Now safe to write
% writecell(T, outFile);
% fprintf('Saved table to %s\n', outFile);

%% 

















