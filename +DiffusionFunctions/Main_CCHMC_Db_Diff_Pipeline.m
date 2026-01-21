clc; clear;



excelFile = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CPIR_Diff_Analysis\Req12702_AB_rawDiff_AD_12092025.xlsx';
% excelFile = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CPIR_Diff_Analysis\Req12702_AB_dcmDiff_AD_12092025.xlsx';
mainDir = '\\rds6.chmccorp.cchmc.org\PulMed-54\CPIR_Images_Database';
WoodsDir = '\\Rds6.cchmc.org\pulmed-35\Woods_CPIR_Images';

% === Load Excel Table ===
T = readcell(excelFile);
% excelFile: full path or file name string
% T: table already read from Excel (e.g., readtable)

[~, excelName, excelExt] = fileparts(excelFile);

if contains(excelName, "rawDiff")
    % Extract relevant columns into simple cell arrays
    StudyCol     = T(:,1);
    SubjectCol   = T(:,2);
    ScanDateCol  = T(:,3);
    SubNumCol    = T(:,4);
    DiffFileCol  = T(:,5);
    ACQ_TypeCol  = T(:,6);
    % SerNum     = T(:,7);
    ScannerSW    = T(:,8);
    ScannerCol   = ''; % T(:,9)   % keep as-is if you truly want blank
    SexCol       = T(:,10);
    AgeCol       = T(:,11);
    DiseaseCol   = T(:,12);
    NoteCol      = T(:,14);
    ImageQCol      = T(:,15);
    RuneCol      = T(:,16);

elseif contains(excelName, "dcmDiff")
    % Extract relevant columns into simple cell arrays
    StudyCol     = T(:,1);
    SubjectCol   = T(:,2);
    ScanDateCol  = T(:,3);
    SubNumCol    = T(:,4);
    DiffFileCol  = T(:,5);
    ACQ_TypeCol  = T(:,6);
    SerNum       = T(:,7);
    ScannerSW    = T(:,8);
    ScannerCol   = T(:,9);
    SexCol       = T(:,12);
    AgeCol       = T(:,13);
    DiseaseCol   = T(:,14);
    NoteCol      = T(:,17);
    RuneCol      = T(:,18);

else
    error("excelFile name must contain 'rawDiff' or 'dcmDiff'. Got: %s", excelFile);
end


nSubjects = size(SexCol,1);

%% 

clc;

for i = 411%:nSubjects % always start from 2
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
    if ismissing(ImageQCol{i})
        ImageQCol{i} = "";
    end
    
    MainInput = struct();
    MainInput.SubjectID        = SubjectCol{i};
    MainInput.Age              = AgeCol{i};
    MainInput.Sex              = SexCol{i}; 
    MainInput.Disease          = DiseaseCol{i}; 
    MainInput.ScanDate         = ScanDateCol{i};
    MainInput.ScannerSoftware  = ScannerSW{i};
    MainInput.SequenceType     = '2D GRE';
    MainInput.denoiseXe        = 'no';
    MainInput.Analyst          = 'Database';
    MainInput.N4Bias           = 'yes';
    MainInput.AnalysisMethod   = 'W.Linear';
    MainInput.AgeCor           = 'no';    
    MainInput.Note             = NoteCol{i};
    MainInput.ImageQuality     = ImageQCol{i};
    MainInput.SliceOrientation = 'transversal';
    MainInput.DiffAcqOrder = 'b-value interleave';
    % diff file
    if contains(DiffFileCol{i}, mainDir, 'IgnoreCase', true) || contains(DiffFileCol{i}, WoodsDir, 'IgnoreCase', true)
        MainInput.diff_file = DiffFileCol{i};
    else
        MainInput.diff_file = fullfile(mainDir, DiffFileCol{i});
    end

    MainInput.diff_file = char(MainInput.diff_file);

    MainInput.ScannerSoftware = ScannerSW{i};
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
    [~, ~, ext] = fileparts(MainInput.diff_file);      
    if strcmpi(ext,'.dcm')
        MainInput.ReconType = 'online';  
        Dinfo = dicominfo(MainInput.diff_file);
        % Image = double(dicomread(MainInput.diff_file));
        StudyTime = Dinfo.StudyTime;
        if strcmpi(ACQ_TypeCol{i},'cpir_diffusion_4bs')
            num_extra_attr1 = 4;
        elseif strcmpi(ACQ_TypeCol{i},'cpir_diffusion_3bs')
            num_extra_attr1 = 3;
        elseif strcmpi(ACQ_TypeCol{i},'cpir_diffusion')
            num_extra_attr1 = 3;
        elseif strcmpi(ACQ_TypeCol{i},'new*_hng_xe_diffusion_29jul2016')
            num_extra_attr1 = 5;
        else

        end
        MainInput.Scanner = 'Philips'; 
        MainInput.sernum  = int2str(SerNum{i});

    elseif strcmp(ext,'.data') 
        MainInput.ReconType = 'offline';
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
        
        diffFilePath = MainInput.diff_file;
        
        % Split path, base name, and extension
        [folderPath, baseName, ~] = fileparts(diffFilePath);
        
        % Construct expected filenames
        sinFile  = fullfile(folderPath, [baseName '.sin']);
        jsonFile = fullfile(folderPath, [baseName '.json']);
        
        % Initialize output
        SeriesTime = [];
        
        % Step 1: Check if .sin file exists
        if isfile(sinFile)
        
            % Step 2: Check if matching .json file exists
            if isfile(jsonFile)
        
                % Step 3: Read and parse JSON
                jsonText = fileread(jsonFile);
                jsonData = jsondecode(jsonText);
        
                % Step 4: Extract SeriesTime if present
                if isfield(jsonData, 'SeriesTime')
                    SeriesTime = jsonData.SeriesTime;
                else
                    warning('JSON file found, but "SeriesTime" field is missing.');
                end
        
            else
                warning('Matching .sin file found, but .json file is missing.');
            end
            MainInput.sernum = SeriesTime(1:6);       
        else
            MainInput.sernum = '000000'; dataset_name; %(5:end);
        end

    else
        MainInput.ReconType = '';
    end

    analysisversion = 'diff_v100';
    MainInput.analysisversion  = analysisversion;
       % Format subject number with leading zeros
    subnum = SubNumCol{i};
    if isnumeric(subnum)
        subnum = num2str(subnum, '%04d');
    end

    MainInput.num_b_values = num_extra_attr1;
    MainInput.Nbvalues  = num_extra_attr1;

    MainInput.EncodingOrder = 'linear'; % confirm this

    % change to analysis folder and then run
    analysisfolder = fullfile(mainDir,StudyCol{i}, 'analysis', analysisversion, ...
        ['sub-', subnum], ['ses-', num2str(MainInput.ScanDate)], ['ser-', (MainInput.sernum)]);
    MainInput.analysisfolder = analysisfolder;

    diff_analysis_folder = fullfile(analysisfolder,'Diffusion_Analysis');
    MainInput.diff_analysis_folder = diff_analysis_folder;
    if ~exist(diff_analysis_folder, 'dir')
        mkdir(diff_analysis_folder);
    end
    MainInput.OutputPath = analysisfolder;
    cd(MainInput.analysisfolder);
    

    % Run pipeline
    DiffusionFunctions.CCHMC_Db_Diff_Pipeline(MainInput)
    % if ismissing(DiffFileCol{i})
    %     VentilationFunctions.CCHMC_Db_Diff_Pipeline_NoData(MainInput)
    % else
    %     if RuneCol{i} == 0
    %         VentilationFunctions.CCHMC_Db_Diff_Pipeline(MainInput);
    %     elseif RuneCol{i} == 1
    %         VentilationFunctions.CCHMC_Db_Diff_Pipeline_rerun(MainInput);
    %     end
    %  
    % end

end