function initial_inputs = read_initial_inputs(excelFile)
%READ_INITIAL_INPUTS Read required batch-analysis inputs from Excel file.

    requiredVars = { ...
        'AnalysisType'
        'SubjectID'
        'SubjectID_PO'
        'Subject'
        'ScanDate'
        'Scan_num'
        'Study'
        'DiseaseType'
        'AGE'
        'SEX'
        'Selected'
        'Recon'
        'ImageQuality'
        'XENON_SERIES_NUMBER'
        'XENON_PATH'
        'PROTON_SERIES_NUMBER'
        'PROTON_PATH'
        'FrequencyOffset_Hz'
        'AnalysisStatus'
        'RunMode'
        'InternalProcessingNotes'
        'PIXEL_SPACING_X'
        'PIXEL_SPACING_Y'
        'SLICE_THICKNESS'
        'Scanner'
        'ScannerSoftware'
        'TRAVERSAL_GEO'
        'SliceOrientation'
        'Polarizer'
        'Note'
        };

    shortVars = { ...
        'Analysis'
        'SubjectID'
        'SubjectID_PO'
        'Subject'
        'Date'
        'Scan'
        'Study'
        'Disease'
        'Age'
        'Sex'
        'Selected'
        'Recon'
        'Quality'
        'XeSeries'
        'XePath'
        'HSeries'
        'HPath'
        'FreqOffset'
        'AnalysisStatus'
        'RunMode'
        'ProcessingNotes'
        'PixelX'
        'PixelY'
        'SliceThickness'
        'Scanner'
        'Software'
        'Traversal'
        'Orientation'
        'Polarizer'
        'Note'
        };

    if ~isfile(excelFile)
        error('Excel file not found:\n%s', excelFile);
    end

    opts = detectImportOptions( ...
        excelFile, ...
        'Sheet', 'all', ...
        'VariableNamingRule', 'preserve');

    data = readtable( ...
        excelFile, ...
        opts, ...
        'Sheet', 'all');

    availableVars = data.Properties.VariableNames;

    missingVars = requiredVars( ...
        ~ismember(requiredVars, availableVars));

    if ~isempty(missingVars)
        fprintf('\nMissing variables:\n');
        fprintf('  %s\n', missingVars{:});

        error( ...
            'The Excel file is missing %d required variable(s).', ...
            numel(missingVars));
    end

    % Extract required variables
    initial_inputs = data(:, requiredVars);

    % Keep only Selected = TRUE
    selected = initial_inputs.Selected;

    if islogical(selected)
        keepRows = selected;
    elseif isnumeric(selected)
        keepRows = selected == 1;
    else
        selected = lower(strtrim(string(selected)));
        keepRows = selected == "true" | selected == "1";
    end

    initial_inputs = initial_inputs(keepRows, :);

    % Rename variables
    initial_inputs.Properties.VariableNames = shortVars;
    
    % ================================================================
    % Add .data extension to offline raw-data paths
    % ================================================================
    for i = 1:height(initial_inputs)
    
        if strcmpi(strtrim(string(initial_inputs.Recon(i))), "Offline")
    
            % Xenon path
            if ~ismissing(initial_inputs.XePath(i)) && ...
                    strlength(strtrim(string(initial_inputs.XePath(i)))) > 0
    
                xePath = char(string(initial_inputs.XePath(i)));
    
                if endsWith(xePath, '.*')
                    xePath = [xePath(1:end-2) '.data'];
                else
                    [~, ~, ext] = fileparts(xePath);
    
                    if isempty(ext)
                        xePath = [xePath '.data'];
                    end
                end
    
                initial_inputs.XePath{i} = xePath;
            end
    
            % Proton path
            if ~ismissing(initial_inputs.HPath(i)) && ...
                    strlength(strtrim(string(initial_inputs.HPath(i)))) > 0
    
                hPath = char(string(initial_inputs.HPath(i)));
    
                if endsWith(hPath, '.*')
                    hPath = [hPath(1:end-2) '.data'];
                else
                    [~, ~, ext] = fileparts(hPath);
    
                    if isempty(ext)
                        hPath = [hPath '.data'];
                    end
                end
    
                initial_inputs.HPath{i} = hPath;
            end
        end
    end
    
    % ================================================================
    % Replace missing Frequency Offset values with 0
    % ================================================================
    if isnumeric(initial_inputs.FreqOffset)
    
        initial_inputs.FreqOffset( ...
            isnan(initial_inputs.FreqOffset)) = 0;
    
    elseif iscell(initial_inputs.FreqOffset)
    
        for i = 1:height(initial_inputs)
    
            value = initial_inputs.FreqOffset{i};
    
            if isempty(value) || ...
                    ismissing(string(value)) || ...
                    strlength(strtrim(string(value))) == 0
    
                initial_inputs.FreqOffset{i} = 0;
            end
        end
    
    else
    
        freqOffset = string(initial_inputs.FreqOffset);
    
        missingFreq = ...
            ismissing(freqOffset) | ...
            strlength(strtrim(freqOffset)) == 0;
    
        freqOffset(missingFreq) = "0";
    
        initial_inputs.FreqOffset = freqOffset;
    end

    % Add row number AFTER filtering
    Row = (1:height(initial_inputs))';

    initial_inputs = addvars( ...
        initial_inputs, ...
        Row, ...
        'Before', 1);

    fprintf('\nInitial inputs loaded successfully.\n');
    fprintf('File          : %s\n', excelFile);
    fprintf('Selected rows : %d\n', height(initial_inputs));
    fprintf('Variables     : %d\n\n', width(initial_inputs));

end