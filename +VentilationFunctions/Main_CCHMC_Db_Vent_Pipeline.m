clc; clear;

%excelFile = '\\rds6.cchmc.org\PulMed-43\CPIR_Share\Carter\08_Master VDP Database Inputs Table\Database_VDP_Inputs_CBM.xlsx';
excelFile = 'C:\Users\MCM5BK\Downloads\VDP_Inputs.xlsx';
%excelFile = 'C:\Users\MCM5BK\OneDrive - cchmc\Documents\03_Data Analysis\02_Data Logs\VDP_Inputs_OfflineCBM.xlsx';

mainDir = '\\Rds6.cchmc.org\PulMed-54\CPIR_Images_Database';
WoodsDir = '\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images';

%%
% === Load Excel Table ===
T = VentilationFunctions.LoadTable(excelFile);

T = VentilationFunctions.enforceCanonicalColumns(T);
VentilationFunctions.validateStandardTable(T);

%% filter T before running the pipeline
T = T(~isnan(cell2mat(T.RunMode)) & cell2mat(T.RunMode) ~= 2, :);
if strcmp(VentilationFunctions.detectAnalyzer(),'CBM')
    T = T(string(T.Polarizer) == "McMaster, Carter", :);
    %T = T(string(T.Polarizer) == "Munidasa, Samal", :);
end
%%
% Extract relevant columns into simple cell arrays
SexCol        = T.Sex;
AgeCol        = T.AGE;
SubjectCol    = T.SubjectID;
DiseaseCol    = T.DiseaseType;
ScanDateCol   = T.ScanDate;
SoftCol       = T.ScannerSoftware;
VoxelX        = T.PIXEL_SPACING_X;
VoxelY        = T.PIXEL_SPACING_Y;
VoxelZ        = T.SLICE_THICKNESS;
StudyCol      = T.Study;
SubNumCol     = T.Subject;
SesNumCol     = T.ScanDate;
ScanNumCol    = num2cell(T.Scan_num);
SeqTypeCol    = T.TRAVERSAL_GEO;
SliceOrient   = T.SliceOrientation;
XeSerNumCol   = T.XeSeriesNumber;
VentFileCol   = T.Vent_path;
HSerNumCol    = T.ProtonSeriesNumber;
AnatFileCol   = T.PROT_FILEPATH;
FrqOffsetCol  = T.FrequencyOffsetHz;
ImageQCol     = T.ImageQuality;
NoteCol       = T.Notes;
RuneCol       = T.RunMode;

nSubjects = size(SexCol,1);

%% 

clc;
for i = 1:nSubjects
    fprintf('Processing subject %d of %d\n', i, nSubjects);

    if ismissing(AgeCol{i})
        AgeCol{i} = "";
    end
    if ismissing(SexCol{i})
        SexCol{i} = "";
    end  
    if ismissing(HSerNumCol{i})
        HSerNumCol{i} = "";
    end          
    if ismissing(AnatFileCol{i})
        AnatFileCol{i} = "";
    end    
    if ismissing(ImageQCol{i})
        ImageQCol{i} = "";
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
    MainInput.Scanner          = 'Philips-3T';
    MainInput.ScannerSoftware  = SoftCol{i};
    MainInput.SequenceType     = SeqTypeCol{i};
    MainInput.denoiseXe        = 'no';
    MainInput.Analyst          = VentilationFunctions.detectAnalyzer();
    MainInput.VoxelSize        = sprintf('[%s,%s,%s]', num2str(VoxelX{i}), num2str(VoxelY{i}), num2str(VoxelZ{i}));
    MainInput.PIXEL_SPACING_X  = VoxelX{i};
    MainInput.PIXEL_SPACING_Y  = VoxelY{i};
    MainInput.SLICE_THICKNESS  = VoxelZ{i};
    MainInput.SliceOrientation = SliceOrient{i};

    % Vent file
    if ismissing(VentFileCol{i})
        MainInput.vent_file = '';
        MainInput.anat_file = '';
        ImageQCol{i} = '1-Failed';
    else
        if contains(fullfile(VentFileCol{i}), mainDir(17:end), 'IgnoreCase', true) || contains(VentFileCol{i}, WoodsDir, 'IgnoreCase', true)
            MainInput.vent_file = fullfile(VentFileCol{i});
        else
            MainInput.vent_file = fullfile(mainDir, VentFileCol{i});
        end
        
        % Anat file
        if contains(fullfile(AnatFileCol{i}), mainDir(17:end), 'IgnoreCase', true) || contains(AnatFileCol{i}, WoodsDir, 'IgnoreCase', true)
            MainInput.anat_file = fullfile(AnatFileCol{i});
        else
            MainInput.anat_file = fullfile(mainDir, AnatFileCol{i});
        end
    end

    MainInput.SCAN_NUM      = ScanNumCol{i};
    MainInput.vent_file     = char(MainInput.vent_file);
    MainInput.anat_file     = char(MainInput.anat_file);

    [~, ~, ext] = fileparts(MainInput.vent_file);      
    if strcmpi(ext,'.dcm')
        MainInput.ReconType = 'online';      
    elseif strcmp(ext,'.data') 
        MainInput.ReconType = 'offline';
    else
        MainInput.ReconType = '';
    end
    analysisversion = 'vent_v100';
    MainInput.analysisversion = analysisversion;

    % Format subject number with leading zeros
    subnum = SubNumCol{i};
    if isnumeric(subnum)
        subnum = num2str(subnum, '%04d');
    end
    
    MainInput.xe_sernum        = XeSerNumCol{i}; %num2str(str2double(OfflineSerNumCol{i}) + 50);
    MainInput.proton_sernum    = HSerNumCol{i};
    MainInput.ImageQuality     = ImageQCol{i};
    MainInput.Note             = NoteCol{i};
    MainInput.freqoffset       = FrqOffsetCol{i};

    MainInput.analysisFolder = fullfile(mainDir, StudyCol{i}, 'analysis', analysisversion, ...
        ['sub-', subnum], ['ses-', num2str(SesNumCol{i})], ['ser-', num2str(MainInput.xe_sernum)]);

    if ~exist(MainInput.analysisFolder, 'dir')
        mkdir(MainInput.analysisFolder);
    end
    cd(MainInput.analysisFolder);

    % Run pipeline
    if ismissing(VentFileCol{i})
        VentilationFunctions.CCHMC_Db_Vent_Pipeline_NoData(MainInput)
    else
        if RuneCol{i} == 0
            VentilationFunctions.CCHMC_Db_Vent_Pipeline(MainInput);
        elseif RuneCol{i} == 1
            VentilationFunctions.CCHMC_Db_Vent_Pipeline_rerun(MainInput);
        elseif RuneCol{i} == 3
            VentilationFunctions.CCHMC_Db_Vent_Pipeline_rerun_registration(MainInput);
        end
        % T{i,24} = MainInput.analysisFolder;
    end
end
% % Define output file path
% [folder, name, ~] = fileparts(excelFile);
% outFile = fullfile(folder, [name '_4.xlsx']);
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

















