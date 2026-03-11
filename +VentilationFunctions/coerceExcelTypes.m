function T = coerceExcelTypes(T)
    import VentilationFunctions.*
    vars = T.Properties.VariableNames;   % cell array of char

    targetStringVars = { ...
        'SubjectID','SubjectID_PO','ScanDate','Session', ...
        'Recon','RunMode','Polarizer','Notes'};

    tf = ismember(vars, targetStringVars);
    stringVars = vars(tf);

    for i = 1:numel(stringVars)
        vn = stringVars{i};
        T.(vn) = string(T.(vn));
        T.(vn) = strtrim(T.(vn));
    end

    targetNumericVars = { ...
        'XeSeriesNumber','ProtonSeriesNumber','FrequencyOffsetHz'};

    tf = ismember(vars, targetNumericVars);
    numericVars = vars(tf);

    for i = 1:numel(numericVars)
        vn = numericVars{i};
        T.(vn) = coerceToDouble(T.(vn));
    end

    if ismember('Selected', vars)
        T.Selected = coerceToLogical(T.Selected);
    end

    if ismember('ScanDate', vars)
        T.ScanDate = normalizeScanDate(T.ScanDate);
    end
end