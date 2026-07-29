clc;
clear;
close all;

main_path = 'C:\Users\bda5ik\Downloads\HC'; % change the main folder to your folder 

correctionFiles = {
    'Fully_sampled.nii.gz',      'FullySampled'
    'SVDKS_correction.nii.gz',   'SVDKS'
    'Global_correction.nii.gz',  'Global'
    'MMCIS_correction.nii.gz',   'MMCIS'
};

maskFileName = 'Mask.nii.gz';
analysisFolderName = 'VDP_analysis';
excelOutputFile = fullfile(main_path, 'All_Subjects_Correction_SNR_VHI_VDP.xlsx');
errorLogFile = fullfile(main_path, 'Processing_Errors.txt');

VentilationTemplate.HeterogeneityIndex = 'yes';
VentilationTemplate.ThreshAnalysis = 'yes';

if exist(excelOutputFile, 'file')
    delete(excelOutputFile);
end

if exist(errorLogFile, 'file')
    delete(errorLogFile);
end

subjectFolders = dir(main_path);
subjectFolders = subjectFolders([subjectFolders.isdir]);
subjectFolders = subjectFolders(~ismember({subjectFolders.name}, {'.', '..'}));

allResults = table();
allErrors = table();

for subjectIndex = 1:numel(subjectFolders)
    subjectID = subjectFolders(subjectIndex).name;
    subjectFolder = fullfile(main_path, subjectID);
    correctionFolder = fullfile(subjectFolder, 'corrections');

    if ~isfolder(correctionFolder)
        fprintf('Skipping %s: corrections folder was not found.\n', subjectID);
        continue;
    end

    maskPath = fullfile(correctionFolder, maskFileName);

    if ~isfile(maskPath)
        fprintf('Skipping %s: %s was not found.\n', subjectID, maskFileName);
        allErrors = [allErrors; table(string(subjectID), "All", "Missing mask file", ...
            'VariableNames', {'SubjectID','CorrectionMethod','ErrorMessage'})];
        continue;
    end

    fprintf('\n============================================================\n');
    fprintf('Processing subject %d of %d: %s\n', subjectIndex, numel(subjectFolders), subjectID);
    fprintf('============================================================\n');

    analysisFolder = fullfile(correctionFolder, analysisFolderName);

    if ~isfolder(analysisFolder)
        mkdir(analysisFolder);
    end

    try
        [lungMask, ~] = loadAndOrientNifti_mask(maskPath);
        lungMask = double(lungMask > 0);
    catch ME
        fprintf('Failed to load mask for %s: %s\n', subjectID, ME.message);
        allErrors = [allErrors; table(string(subjectID), "All", string(ME.message), ...
            'VariableNames', {'SubjectID','CorrectionMethod','ErrorMessage'})];
        continue;
    end

    Outputs = struct();
    Outputs.SubjectID = subjectID;
    Outputs.CorrectionFolder = correctionFolder;
    Outputs.AnalysisFolder = analysisFolder;
    Outputs.MaskFile = maskPath;
    Outputs.LungMask = lungMask;
    Outputs.ProcessingDate = datestr(now, 30);

    subjectResults = table();


    for methodIndex = 1:size(correctionFiles, 1)
        correctionFileName = correctionFiles{methodIndex, 1};
        methodName = correctionFiles{methodIndex, 2};
        correctionPath = fullfile(correctionFolder, correctionFileName);
        methodField = matlab.lang.makeValidName(methodName);

        fprintf('\n%s | %s\n', subjectID, correctionFileName);

        if ~isfile(correctionPath)
            warning('%s was not found for subject %s.', correctionFileName, subjectID);
            allErrors = [allErrors; table(string(subjectID), string(methodName), ...
                "Missing correction image", ...
                'VariableNames', {'SubjectID','CorrectionMethod','ErrorMessage'})];
            continue;
        end

        try
            [MR, imageReferenceInfo] = loadAndOrientNifti(correctionPath);
            MR = double(MR);

            if ~isequal(size(MR), size(lungMask))
                error('Image size %s does not match mask size %s.', ...
                    mat2str(size(MR)), mat2str(size(lungMask)));
            end

            Ventilation = VentilationTemplate;
            Ventilation.parentPath = analysisFolder;
            Ventilation.Image = MR;
            Ventilation.UncorrectedImage = MR;
            Ventilation.Image_uncorrected = MR;
            Ventilation.LungMask = lungMask;
            Ventilation.AirwayMask = zeros(size(lungMask), 'double');
            Ventilation.VesselMask = zeros(size(lungMask), 'double');

            MainInput = struct();
            MainInput.XeFileName = correctionFileName;
            MainInput.XeDataLocation = correctionFolder;
            MainInput.XeFullPath = correctionPath;

            Proton = struct();
            Proton.Image = zeros(size(lungMask), 'double');
            Proton.ProtonRegistered = zeros(size(lungMask), 'double');
            Proton.ProtonRegisteredColored = zeros(size(lungMask), 'double');
            methodOutputFolder = fullfile(analysisFolder, methodName);

            if ~isfolder(methodOutputFolder)
                mkdir(methodOutputFolder);
            end

            Outputs.(methodField) = struct();
            Outputs.(methodField).InputFile = correctionPath;
            Outputs.(methodField).ImageReferenceInfo = imageReferenceInfo;
            Outputs.(methodField).Uncorrected = struct();
            Outputs.(methodField).N4 = struct();

            Ventilation.Image = MR;
            Ventilation.UncorrectedImage = MR;
            Ventilation.LungMask = lungMask;
            Ventilation.VesselMask = zeros(size(lungMask), 'double');

            [Ventilation] = VentilationFunctions.calculate_SNR(Ventilation);
            close all;

            Outputs.(methodField).SNR_slice = Ventilation.SNR_slice;
            Outputs.(methodField).SNRvv_slice = Ventilation.SNRvv_slice;
            Outputs.(methodField).SNR_lung = Ventilation.SNR_lung;
            Outputs.(methodField).SNR_vv = Ventilation.SNR_vv;

            Outputs.(methodField).Uncorrected.Image = MR;

            if strcmpi(Ventilation.HeterogeneityIndex, 'yes')
                [Ventilation] = VentilationFunctions.calculate_VHI(Ventilation);
                close all;

                Outputs.(methodField).Uncorrected.CV_maps = Ventilation.CV_maps;
                Outputs.(methodField).Uncorrected.sliceMeanCV = Ventilation.sliceMeanCV;
                Outputs.(methodField).Uncorrected.sliceVHI = Ventilation.sliceVHI;
                Outputs.(methodField).Uncorrected.overallMeanCV = Ventilation.overallMeanCV;
                Outputs.(methodField).Uncorrected.overallVHI = Ventilation.overallVHI;

                saveMapIfPossible(Ventilation.CV_maps, ...
                    fullfile(methodOutputFolder, [methodName '_Uncorrected_CV_Map.nii.gz']), ...
                    imageReferenceInfo);
            end

            maskarray = double(Ventilation.LungMask + Ventilation.VesselMask) .* Ventilation.LungMask;
            maskarray(maskarray > 1) = 0;
            maskarray = double(maskarray);
            Ventilation.LungMask = maskarray;
            Ventilation.N4Analysis = 'no';
            Ventilation.IncompleteThresh = 60;
            Ventilation.CompleteThresh = 30;
            Ventilation.HyperventilatedThresh = 200;
            Ventilation.MedianFilter = 'yes';
            MainInput.NoProtonImage = 'yes';
            if strcmpi(Ventilation.ThreshAnalysis, 'yes')
                [Ventilation] = VentilationFunctions.calculate_VDP_CCHMC(Ventilation, Proton, MainInput);
                close all;

                Outputs.(methodField).Uncorrected.TH60.VDP = Ventilation.Threshold.VDP;
                Outputs.(methodField).Uncorrected.TH60.defectArray = Ventilation.Threshold.defectArray;
                Outputs.(methodField).Uncorrected.TH60.BinPrecentValues = Ventilation.Threshold.THBins;

                saveMapIfPossible(Ventilation.Threshold.defectArray, ...
                    fullfile(methodOutputFolder, [methodName '_Uncorrected_TH60_Defect_Map.nii.gz']), ...
                    imageReferenceInfo);
            end

            saveMapIfPossible(MR, ...
                fullfile(methodOutputFolder, [methodName '_Uncorrected_Image.nii.gz']), ...
                imageReferenceInfo);

            n4WorkingFolder = fullfile(methodOutputFolder, 'N4_working');

            if ~isfolder(n4WorkingFolder)
                mkdir(n4WorkingFolder);
            end

            [N4, BiasField] = VentilationFunctions.N4_bias_correction(MR, maskarray, n4WorkingFolder);
            close all;

            MR_N4 = double(N4);

            Ventilation.Image_uncorrected = MR;
            Ventilation.Image = MR_N4;
            Ventilation.UncorrectedImage = MR;
            Ventilation.LungMask = maskarray;
            Ventilation.VesselMask = zeros(size(maskarray), 'double');
            Ventilation.parentPath = n4WorkingFolder;

            Outputs.(methodField).N4.Image = MR_N4;
            Outputs.(methodField).N4.BiasField = BiasField;

            if strcmpi(Ventilation.HeterogeneityIndex, 'yes')
                [Ventilation] = VentilationFunctions.calculate_VHI(Ventilation);
                close all;

                Outputs.(methodField).N4.CV_maps = Ventilation.CV_maps;
                Outputs.(methodField).N4.sliceMeanCV = Ventilation.sliceMeanCV;
                Outputs.(methodField).N4.sliceVHI = Ventilation.sliceVHI;
                Outputs.(methodField).N4.overallMeanCV = Ventilation.overallMeanCV;
                Outputs.(methodField).N4.overallVHI = Ventilation.overallVHI;

                saveMapIfPossible(Ventilation.CV_maps, ...
                    fullfile(methodOutputFolder, [methodName '_N4_CV_Map.nii.gz']), ...
                    imageReferenceInfo);
            end

            if strcmpi(Ventilation.ThreshAnalysis, 'yes')
                [Ventilation] = VentilationFunctions.calculate_VDP_CCHMC(Ventilation, Proton, MainInput);
                close all;

                Outputs.(methodField).N4.TH60.VDP = Ventilation.Threshold.VDP;
                Outputs.(methodField).N4.TH60.defectArray = Ventilation.Threshold.defectArray;
                Outputs.(methodField).N4.TH60.BinPrecentValues = Ventilation.Threshold.THBins;

                saveMapIfPossible(Ventilation.Threshold.defectArray, ...
                    fullfile(methodOutputFolder, [methodName '_N4_TH60_Defect_Map.nii.gz']), ...
                    imageReferenceInfo);
            end

            saveMapIfPossible(MR_N4, ...
                fullfile(methodOutputFolder, [methodName '_N4_Image.nii.gz']), ...
                imageReferenceInfo);

            if isnumeric(BiasField) && isequal(size(BiasField), size(MR))
                saveMapIfPossible(BiasField, ...
                    fullfile(methodOutputFolder, [methodName '_N4_Bias_Field.nii.gz']), ...
                    imageReferenceInfo);
            end

            resultRow = table();
            resultRow.SubjectID = string(subjectID);
            resultRow.CorrectionMethod = string(methodName);
            resultRow.InputFile = string(correctionFileName);
            resultRow.SNR_Lung = scalarOrNaN(Outputs.(methodField).SNR_lung);
            resultRow.SNR_VentilatedVolume = scalarOrNaN(Outputs.(methodField).SNR_vv);
            resultRow.Uncorrected_OverallMeanCV = scalarOrNaN(getNestedField(Outputs.(methodField), {'Uncorrected','overallMeanCV'}));
            resultRow.Uncorrected_OverallVHI = scalarOrNaN(getNestedField(Outputs.(methodField), {'Uncorrected','overallVHI'}));
            resultRow.Uncorrected_TH60_VDP = scalarOrNaN(getNestedField(Outputs.(methodField), {'Uncorrected','TH60','VDP'}));
            resultRow.N4_OverallMeanCV = scalarOrNaN(getNestedField(Outputs.(methodField), {'N4','overallMeanCV'}));
            resultRow.N4_OverallVHI = scalarOrNaN(getNestedField(Outputs.(methodField), {'N4','overallVHI'}));
            resultRow.N4_TH60_VDP = scalarOrNaN(getNestedField(Outputs.(methodField), {'N4','TH60','VDP'}));

            uncorrectedBins = rowVectorOrEmpty(getNestedField(Outputs.(methodField), {'Uncorrected','TH60','BinPrecentValues'}));
            n4Bins = rowVectorOrEmpty(getNestedField(Outputs.(methodField), {'N4','TH60','BinPrecentValues'}));

            for binIndex = 1:max(numel(uncorrectedBins), numel(n4Bins))
                uncorrectedVariable = sprintf('Uncorrected_TH60_Bin_%d_Percent', binIndex);
                n4Variable = sprintf('N4_TH60_Bin_%d_Percent', binIndex);

                if binIndex <= numel(uncorrectedBins)
                    resultRow.(uncorrectedVariable) = uncorrectedBins(binIndex);
                else
                    resultRow.(uncorrectedVariable) = NaN;
                end

                if binIndex <= numel(n4Bins)
                    resultRow.(n4Variable) = n4Bins(binIndex);
                else
                    resultRow.(n4Variable) = NaN;
                end
            end

            subjectResults = appendCompatibleTables(subjectResults, resultRow);
            allResults = appendCompatibleTables(allResults, resultRow);

            methodMatFile = fullfile(methodOutputFolder, [methodName '_Results.mat']);
            MethodOutputs = Outputs.(methodField);
            save(methodMatFile, 'MethodOutputs', '-v7.3');

            fprintf('Completed %s for %s.\n', methodName, subjectID);

        catch ME
            close all;
            warning('Failed %s for %s: %s', methodName, subjectID, ME.message);

            allErrors = [allErrors; table(string(subjectID), string(methodName), string(ME.message), ...
                'VariableNames', {'SubjectID','CorrectionMethod','ErrorMessage'})];

            Outputs.(methodField).Error = struct();
            Outputs.(methodField).Error.Message = ME.message;
            Outputs.(methodField).Error.Identifier = ME.identifier;
            Outputs.(methodField).Error.Stack = ME.stack;
        end
    end

    outputsMatFile = fullfile(subjectFolder, 'Outputs.mat');
    save(outputsMatFile, 'Outputs', '-v7.3');

    if ~isempty(subjectResults)
        subjectExcelFile = fullfile(subjectFolder, [subjectID '_Correction_SNR_VHI_VDP.xlsx']);

        if isfile(subjectExcelFile)
            delete(subjectExcelFile);
        end

        writetable(subjectResults, subjectExcelFile, 'Sheet', 'Results');
    end

    fprintf('\nSaved subject outputs:\n%s\n', outputsMatFile);
end

if ~isempty(allResults)
    writetable(allResults, excelOutputFile, 'Sheet', 'All Results');

    methodSummary = groupsummary(allResults, 'CorrectionMethod', {'mean','std'}, ...
        {'SNR_Lung','SNR_VentilatedVolume', ...
         'Uncorrected_OverallMeanCV','Uncorrected_OverallVHI','Uncorrected_TH60_VDP', ...
         'N4_OverallMeanCV','N4_OverallVHI','N4_TH60_VDP'});

    writetable(methodSummary, excelOutputFile, 'Sheet', 'Method Summary');
end

if ~isempty(allErrors)
    if ~isfile(excelOutputFile)
        writetable(allErrors, excelOutputFile, 'Sheet', 'Errors');
    else
        writetable(allErrors, excelOutputFile, 'Sheet', 'Errors');
    end

    fileID = fopen(errorLogFile, 'w');

    if fileID ~= -1
        for errorIndex = 1:height(allErrors)
            fprintf(fileID, '%s | %s | %s\n', ...
                allErrors.SubjectID(errorIndex), ...
                allErrors.CorrectionMethod(errorIndex), ...
                allErrors.ErrorMessage(errorIndex));
        end

        fclose(fileID);
    end
end

fprintf('\n============================================================\n');
fprintf('Processing complete.\n');
fprintf('Combined Excel file:\n%s\n', excelOutputFile);
fprintf('============================================================\n');

function [volume, referenceInfo] = loadAndOrientNifti(filePath)
    try
        loadedNifti = LoadData.load_nii(filePath);
    catch
        loadedNifti = LoadData.load_untouch_nii(filePath);
    end

    volume = double(squeeze(loadedNifti.img));
    volume = permute(volume, [1,3,2]);
    volume = imrotate(volume, -90);
    volume = flip(volume, 1);
    imslice(volume)

    try
        referenceInfo = niftiinfo(filePath);
    catch
        referenceInfo = [];
    end
end

function [volume, referenceInfo] = loadAndOrientNifti_mask(filePath)
    try
        loadedNifti = LoadData.load_nii(filePath);
    catch
        loadedNifti = LoadData.load_untouch_nii(filePath);
    end

    volume = double(squeeze(loadedNifti.img));
    volume = permute(volume, [1,3,2]);
    volume = imrotate(volume, -90);
    volume = flip(volume, 1);
    % imslice(volume)

    try
        referenceInfo = niftiinfo(filePath);
    catch
        referenceInfo = [];
    end
end

function saveMapIfPossible(map, outputPath, referenceInfo)
    if isempty(map) || ~isnumeric(map)
        return;
    end

    map = double(squeeze(map));

    if ndims(map) ~= 3
        warning('Map was not saved as NIfTI because it is not 3-D: %s', outputPath);
        return;
    end

    mapOriginalOrientation = imrotate(flip(map, 2), -90);

    [outputFolder, outputName, outputExtension] = fileparts(outputPath);

    if strcmpi(outputExtension, '.gz')
        [~, outputNameWithoutNii, ~] = fileparts(outputName);
        baseOutputPath = fullfile(outputFolder, outputNameWithoutNii);
    else
        baseOutputPath = fullfile(outputFolder, outputName);
    end

    if isempty(referenceInfo)
        niftiwrite(single(mapOriginalOrientation), baseOutputPath, 'Compressed', true);
        return;
    end

    outputInfo = referenceInfo;
    outputInfo.Datatype = 'single';
    outputInfo.ImageSize = size(mapOriginalOrientation);

    try
        niftiwrite(single(mapOriginalOrientation), baseOutputPath, outputInfo, 'Compressed', true);
    catch
        niftiwrite(single(mapOriginalOrientation), baseOutputPath, 'Compressed', true);
    end
end

function value = getNestedField(structure, fieldPath)
    value = [];

    try
        currentValue = structure;

        for pathIndex = 1:numel(fieldPath)
            fieldName = fieldPath{pathIndex};

            if ~isstruct(currentValue) || ~isfield(currentValue, fieldName)
                value = [];
                return;
            end

            currentValue = currentValue.(fieldName);
        end

        value = currentValue;
    catch
        value = [];
    end
end

function value = scalarOrNaN(inputValue)
    if isempty(inputValue) || ~isnumeric(inputValue)
        value = NaN;
        return;
    end

    finiteValues = double(inputValue(isfinite(inputValue)));

    if isempty(finiteValues)
        value = NaN;
    elseif isscalar(finiteValues)
        value = finiteValues;
    else
        value = mean(finiteValues(:), 'omitnan');
    end
end

function outputVector = rowVectorOrEmpty(inputValue)
    if isempty(inputValue) || ~isnumeric(inputValue)
        outputVector = [];
    else
        outputVector = double(inputValue(:).');
    end
end

function outputTable = appendCompatibleTables(existingTable, newRow)
    if isempty(existingTable)
        outputTable = newRow;
        return;
    end

    existingNames = existingTable.Properties.VariableNames;
    newNames = newRow.Properties.VariableNames;
    allNames = union(existingNames, newNames, 'stable');

    for nameIndex = 1:numel(allNames)
        variableName = allNames{nameIndex};

        if ~ismember(variableName, existingNames)
            existingTable.(variableName) = NaN(height(existingTable), 1);
        end

        if ~ismember(variableName, newNames)
            if isstring(existingTable.(variableName))
                newRow.(variableName) = strings(1, 1);
            elseif iscell(existingTable.(variableName))
                newRow.(variableName) = {[]};
            else
                newRow.(variableName) = NaN;
            end
        end
    end

    existingTable = existingTable(:, allNames);
    newRow = newRow(:, allNames);
    outputTable = [existingTable; newRow];
end
