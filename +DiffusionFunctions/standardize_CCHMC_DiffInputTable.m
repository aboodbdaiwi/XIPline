function C = standardize_CCHMC_DiffInputTable(excelFile)
%STANDARDIZEDIFFINPUTTABLE Read diffusion input Excel file and return
%legacy column arrays expected by Main_CCHMC_Db_Diff_Pipeline.

    raw = readcell(excelFile);

    headers = string(raw(1,:));
    data = raw(2:end,:);   % row 2 is intentionally skipped

    % Decide format from headers
    if any(headers == "AnalysisType") && any(headers == "IMAGE_PATH")
        C = fromWebAppDiffTable(headers, data);
    elseif any(headers == "diff_data_path") || any(headers == "STUDY")
        C = fromRawDiffTable(headers, data);
    else
        error("standardizeDiffInputTable:UnknownFormat", ...
            "Could not identify diffusion input table format: %s", excelFile);
    end

    function C = fromWebAppDiffTable(raw)
        data = raw(2:end,:);   % use raw(3:end,:) only if row 2 is not real data

        % Keep only diffusion rows
        isDiff = strcmpi(string(data(:,1)), "diff");
        data = data(isDiff,:);
    
        n = size(data,1);
    
        C = struct();
    
        C.SubjectCol  = data(:,2);    % SubjectID
        C.SubNumCol   = data(:,4);    % Subject
        C.ScanDateCol = normalizeDateCol(data(:,5));
        C.DiseaseCol  = data(:,8);    % DiseaseType
        C.AgeCol      = data(:,9);    % AGE
        C.SexCol      = normalizeSexCol(data(:,10));
        C.ImageQCol   = data(:,13);   % ImageQuality
        C.SerNum      = data(:,14);   % IMAGE_SERIES_NUMBER
        C.DiffFileCol = normalizeDiffPathCol(data(:,15));
        C.RuneCol     = data(:,22);   % RunMode
        C.ScannerSW   = data(:,28);   % ScannerSoftware
        C.NoteCol     = data(:,32);   % Note
    
        % Derived / default fields expected by legacy diffusion script
        C.StudyCol    = inferStudyCol(C.SubjectCol, C.DiffFileCol);
        C.ACQ_TypeCol = repmat({'cpir_diffusion'}, n, 1);
        C.ScannerCol  = repmat({''}, n, 1);
        C.MaskCol     = repmat({''}, n, 1);
    end

    function C = fromRawDiffTable(headers, data)
    
        C = struct();
    
        C.StudyCol     = getCol(headers, data, ["STUDY", "Study"]);
        C.SubjectCol   = getCol(headers, data, ["SUBJECT_ID", "SubjectID", "Subject ID"]);
        C.ScanDateCol  = getCol(headers, data, ["SCAN_DATE", "ScanDate", "Scan Date"]);
        C.SubNumCol    = getCol(headers, data, ["subject number", "Subject", "SUBJECT"]);
        C.DiffFileCol  = getCol(headers, data, ["diff_data_path", "DiffFile", "DIFF_FILEPATH_NEW"]);
        C.ACQ_TypeCol  = getCol(headers, data, ["ACQ_DESC_LIST", "ACQ_Type", "AcqType"]);
    
        if hasCol(headers, ["SerNum", "SER_NUM", "IMAGE_SERIES_NUMBER"])
            C.SerNum = getCol(headers, data, ["SerNum", "SER_NUM", "IMAGE_SERIES_NUMBER"]);
        else
            C.SerNum = repmat({NaN}, size(data,1), 1);
        end
    
        C.ScannerSW    = getCol(headers, data, ["SW_RELEASE", "ScannerSoftware", "Scanner Software"]);
        C.ScannerCol   = repmat({''}, size(data,1), 1);
    
        C.SexCol       = getCol(headers, data, ["SEX", "Sex"]);
        C.AgeCol       = getCol(headers, data, ["Age", "AGE"]);
        C.DiseaseCol   = getCol(headers, data, ["disease", "Disease", "DiseaseType"]);
        C.NoteCol      = getCol(headers, data, ["Note", "Notes"]);
        C.ImageQCol    = getCol(headers, data, ["IQ", "ImageQuality", "Image Quality"]);
        C.RuneCol      = getCol(headers, data, ["Run", "RunMode"]);
    
        if hasCol(headers, ["Mask_path", "MaskPath", "Mask Path"])
            C.MaskCol = getCol(headers, data, ["Mask_path", "MaskPath", "Mask Path"]);
        else
            C.MaskCol = repmat({''}, size(data,1), 1);
        end
    end
    function col = getCol(headers, data, candidates)

        candidates = string(candidates);
    
        for k = 1:numel(candidates)
            idx = find(headers == candidates(k), 1);
            if ~isempty(idx)
                col = data(:,idx);
                return
            end
        end
    
        error("standardize_CCHMC_DiffInputTable:MissingColumn", ...
            "Missing required column. Tried: %s", strjoin(candidates, ", "));
    end
    
    
    function tf = hasCol(headers, candidates)
    
        tf = false;
    
        for k = 1:numel(candidates)
            if any(headers == candidates(k))
                tf = true;
                return
            end
        end
    end
end