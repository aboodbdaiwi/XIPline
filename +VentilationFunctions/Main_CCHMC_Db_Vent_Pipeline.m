clc; clear;


excelFile = 'C:\Users\MCM5BK\OneDrive - cchmc\Documents\42_VDP Database Updating Project\04_VDP_inputs_log\Database_VDP_Inputs - Copy.xlsx';
mainDir = '\\rds6.chmccorp.cchmc.org\PulMed-54\CPIR_Images_Database';
WoodsDir = '\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images';

% === Load Excel Table ===
T = readcell(excelFile);

% Extract relevant columns into simple cell arrays
SexCol       = T(:,1);
AgeCol       = T(:,2);
SubjectCol   = T(:,3);
DiseaseCol   = T(:,4);
ScanDateCol  = T(:,5);
SoftCol      = T(:,6);
VoxelX       = T(:,7);
VoxelY       = T(:,8);
VoxelZ       = T(:,9);
StudyCol      = T(:,10);
SubNumCol    = T(:,11);
SesNumCol    = T(:,12);
ScanNumCol   = T(:,13);
SeqTypeCol   = T(:,14);
SliceOrient  = T(:,15);
XeSerNumCol    = T(:,16);
VentFileCol  = T(:,17);
HSerNumCol    = T(:,18);
AnatFileCol  = T(:,19);
FrqOffsetCol    = T(:,20);
ImageQCol    = T(:,21);
NoteCol      = T(:,22);
RuneCol      = T(:,23);

nSubjects = size(SexCol,1);

%% 

clc;
for i = 33:nSubjects % always start from 2
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
    MainInput.Analyst          = 'Database';
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
        if contains(VentFileCol{i}, mainDir, 'IgnoreCase', true) || contains(VentFileCol{i}, WoodsDir, 'IgnoreCase', true)
            MainInput.vent_file = VentFileCol{i};
        else
            MainInput.vent_file = fullfile(mainDir, VentFileCol{i});
        end
        
        % Anat file
        if contains(AnatFileCol{i}, mainDir, 'IgnoreCase', true) || contains(AnatFileCol{i}, WoodsDir, 'IgnoreCase', true)
            MainInput.anat_file = AnatFileCol{i};
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

















