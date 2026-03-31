function T = enforceCanonicalColumns(T)
%ENFORCECANONICALCOLUMNS Guarantee consistent variables, order, and presence.

    canonicalVars = [
        "Selected"
        "SubjectID"
        "SubjectID_PO"
        "ScanDate"
        "Session"
        "XeSeriesNumber"
        "ProtonSeriesNumber"
        "FrequencyOffsetHz"
        "Recon"
        "RunMode"
        "Polarizer"
        "Notes"
        "Sex"
    ];

    % Ensure every canonical column exists
    for i = 1:numel(canonicalVars)
        varName = canonicalVars(i);

        if ~ismember(varName, T.Properties.VariableNames)
            switch varName
                case "Selected"
                    T.(varName) = false(height(T),1);

                case {"XeSeriesNumber","ProtonSeriesNumber","FrequencyOffsetHz"}
                    T.(varName) = nan(height(T),1);

                otherwise
                    T.(varName) = strings(height(T),1);
            end
        end
    end

    % Reorder and optionally drop extra columns
    %T = T(:, canonicalVars);
end