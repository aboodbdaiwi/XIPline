clc; clear;


excelFile = 'C:\Users\MCM5BK\OneDrive - cchmc\Documents\56_diffusion pipeline\diff_main_input.xlsx';
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
DiffFileCol   = T(:,10);
NoteCol      = T(:,11);
RuneCol      = T(:,12);

nSubjects = size(SexCol,1);

%% 

clc;

for i = 2:6 % always start from 2
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
    MainInput.SequenceType     = '2D GRE';
    MainInput.denoiseXe        = 'no';
    MainInput.Analyst          = 'Database';
    MainInput.N4Bias           = 'yes';
    MainInput.AnalysisMethod   = '';
    MainInput.AgeCor           = 'no';   
    MainInput.Note             = NoteCol{i};
    MainInput.SliceOrientation = 'transversal';

    % diff file
    if contains(DiffFileCol{i}, mainDir, 'IgnoreCase', true) || contains(DiffFileCol{i}, WoodsDir, 'IgnoreCase', true)
        MainInput.diff_file = DiffFileCol{i};
    else
        MainInput.diff_file = fullfile(mainDir, DiffFileCol{i});
    end

    MainInput.diff_file = char(MainInput.diff_file);

    [~, ~, ext] = fileparts(MainInput.diff_file);      
    if strcmpi(ext,'.dcm')
        MainInput.ReconType = 'online';      
    elseif strcmp(ext,'.data') 
        MainInput.ReconType = 'offline';
    else
        MainInput.ReconType = '';
    end
    analysisversion = 'diff_v100';

       % Format subject number with leading zeros
    subnum = SubNumCol{i};
    if isnumeric(subnum)
        subnum = num2str(subnum, '%04d');
    end

    % extract bvalues and file name from .list file
    toks = regexp(MainInput.diff_file,'^(.*?)(\.list|\.data)?$','tokens');
    prefix = toks{1}{1};
    listname = sprintf('%s.list',prefix);
 
    fid = fopen(listname,'r');

    dataset_name = '';
    num_extra_attr1 = NaN;
    
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        
        % Check for dataset name
        if contains(line, 'Dataset name:')
            tokens = regexp(line, 'Dataset name:\s*(\S+)', 'tokens');
            if ~isempty(tokens)
                dataset_name = tokens{1}{1};
            end
        end
        
        % Check for number_of_extra_attribute_1_values
        if contains(line, 'number_of_extra_attribute_1_values')
            tokens = regexp(line, 'number_of_extra_attribute_1_values\s*:\s*(\d+)', 'tokens');
            if ~isempty(tokens)
                num_extra_attr1 = str2double(tokens{1}{1});
            end
        end
    end
    
    fclose(fid);
    
    MainInput.sernum = dataset_name(5:end);
    MainInput.num_b_values = num_extra_attr1;
 

    MainInput.EncodingOrder = 'linear'; % confirm this

    % change to analysis folder and then run
    analysisfolder = fullfile(mainDir,StudyCol{i}, 'analysis', analysisversion, ...
        ['sub-', subnum], ['ses-', num2str(MainInput.ScanDate)], ['ser-', MainInput.sernum]);
    MainInput.analysisfolder = analysisfolder;

    diff_analysis_folder = fullfile(analysisfolder,'Diffusion_Analysis');
    MainInput.diff_analysis_folder = diff_analysis_folder;
    if ~exist(diff_analysis_folder, 'dir')
        mkdir(diff_analysis_folder);
    end
    cd(MainInput.analysisfolder);

     
    
    
    
    % Run pipeline
    DiffusionFunctions.CCHMC_Db_Diff_Pipeline(MainInput)
    % if ismissing(DiffFileCol{i})
    %     VentilationFunctions.CCHMC_Db_Vent_Pipeline_NoData(MainInput)
    % else
    %     if RuneCol{i} == 0
    %         VentilationFunctions.CCHMC_Db_Vent_Pipeline(MainInput);
    %     elseif RuneCol{i} == 1
    %         VentilationFunctions.CCHMC_Db_Vent_Pipeline_rerun(MainInput);
    %     end
    %  
    % end

end