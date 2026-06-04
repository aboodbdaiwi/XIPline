clc; clear;

% Input
excelFile = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CPIR_Diff_Analysis\Validation\main_rawDiff.xlsx';
sheetNum  = 3;

outputDir = fileparts(excelFile);
outExcel  = fullfile(outputDir,'healthy_reference.xlsx');
outMat    = fullfile(outputDir,'healthy_reference.mat');

% Load Excel
T = readcell(excelFile,'Sheet',sheetNum);

% Column 21 contains analysis folders
analysisFolders = T(2:end,21);

% Initialize output table
healthy_reference = table();

%% Loop through folders
for i = 1:numel(analysisFolders)

    if isempty(analysisFolders{i}) || ~ischar(analysisFolders{i})
        warning('Row %d skipped: empty or invalid folder path.', i+1);
        continue
    end

    ADC_path = fullfile(analysisFolders{i}, ...
        'Diffusion_Analysis','workspace.mat');

    if ~isfile(ADC_path)
        warning('File not found: %s', ADC_path);
        continue
    end

    fprintf('Loading %s\n', ADC_path);

    ADC_sub = load(ADC_path);

    % Basic subject info
    row = table();

    row.SubjectID = string(ADC_sub.MainInput.SubjectID);
    row.Age       = ADC_sub.MainInput.Age;
    row.Sex       = string(ADC_sub.MainInput.Sex);
    row.ScanDate  = string(ADC_sub.MainInput.ScanDate);

    % ADC metrics
    row.meanADC      = ADC_sub.Diffusion.meanADC;
    row.stdADC       = ADC_sub.Diffusion.stdADC;
    row.ADC_cv       = ADC_sub.Diffusion.stdADC / ADC_sub.Diffusion.meanADC;
    row.ADC_skewness = ADC_sub.Diffusion.ADC_skewness;
    row.ADC_kurtosis = ADC_sub.Diffusion.ADC_kurtosis;

    % R2 metrics excluding zeros
    if isfield(ADC_sub.Diffusion,'Rmap')
        R2vals = ADC_sub.Diffusion.Rmap(ADC_sub.Diffusion.Rmap ~= 0 & ~isnan(ADC_sub.Diffusion.Rmap));

        if ~isempty(R2vals)
            row.mean_R2 = mean(R2vals);
            row.sd_R2   = std(R2vals);
        else
            row.mean_R2 = NaN;
            row.sd_R2   = NaN;
        end
    else
        row.mean_R2 = NaN;
        row.sd_R2   = NaN;
    end

    A = ADC_sub.Diffusion.LB_colorDiffmap;
    I = zeros(size(A,1), size(A,2), size(A,3));
    for s = 1:size(A,3)
        RGBslice = squeeze(A(:,:,s,:));
        I(:,:,s) = rgb2gray(RGBslice);
    end
    ADC_sub.Diffusion.LB_ADCmap = I;
    % Morphometry metrics
    morphFields = { ...
        'R_mean','R_mean','h_mean','r_mean','Lm_mean','SVR_mean','Na_mean', ...
        'R_std','h_std','r_std','Lm_std','SVR_std','Na_std', ...
        'DDC_mean','alpha_mean','LmD_mean', ...
        'DDC_std','alpha_std','LmD_std','LDR','NDR','HDR'};

    for f = 1:numel(morphFields)
        fieldName = morphFields{f};

        if isfield(ADC_sub.Diffusion,fieldName)
            row.(fieldName) = ADC_sub.Diffusion.(fieldName);
        else
            row.(fieldName) = NaN;
        end
    end

    % Store maps separately, not in Excel
    mapFields = { ...
        'ADCmap','R_map','h_map','r_map','Lm_map','SVR_map','Na_map','So_map', ...
        'DDC_map','alpha_map','SEMSo_map','LmD_map','LB_ADCmap'};

    maps(i).SubjectID = row.SubjectID; %#ok<SAGROW>
    maps(i).ADC_path  = ADC_path;

    for f = 1:numel(mapFields)
        fieldName = mapFields{f};

        if isfield(ADC_sub.Diffusion,fieldName)
            maps(i).(fieldName) = ADC_sub.Diffusion.(fieldName);
        else
            maps(i).(fieldName) = [];
        end
    end

    % Append row
    healthy_reference = [healthy_reference; row];

end

% Save outputs
writetable(healthy_reference,outExcel);
save(outMat,'healthy_reference','maps');

fprintf('\nSaved:\n%s\n%s\n',outExcel,outMat);