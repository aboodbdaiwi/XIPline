function validateStandardTable(T)
%VALIDATEVDPTABLE Fail early if required canonical fields are missing or bad.

    requiredVars = [
        "SubjectID"
        "ScanDate"
        "ProtonSeriesNumber"
    ];

    for i = 1:numel(requiredVars)
        if ~ismember(requiredVars(i), T.Properties.VariableNames)
            error("validateVDPTable:MissingColumn", ...
                "Missing required column: %s", requiredVars(i));
        end
    end

    % Example sanity checks
    if ~islogical(T.Selected)
        error("validateVDPTable:BadType", ...
            "Selected must be logical after standardization.");
    end

    if ~(iscell(T.ProtonSeriesNumber) && all(cellfun(@isnumeric, T.ProtonSeriesNumber)))
        error("validateVDPTable:BadType", ...
            "ProtonSeriesNumber must be numeric after standardization.");
    end
end